function synthetic_kinetic_experiment(;
    identifier="synthetic_kinetic_rate10",
    heating_rate_K_per_min=10.0,
    activation_energy_kJ_per_mol=145.0,
    preexponential_per_min=1.0e9,
    m=0.25,
    n=1.10,
    p=0.0,
    model=:sestak_berggren_2,
    noisy=false,
    inject_nonpositive=false,
    composition=Composition(0.5),
)
    temperature_K = collect(450.0:0.5:950.0)
    alpha = zeros(Float64, length(temperature_K))
    alpha[1] = 1.0e-3
    beta = Float64(heating_rate_K_per_min)
    kinetic_rate(temperature, conversion) =
        preexponential_per_min *
        exp(
            -activation_energy_kJ_per_mol /
            (IsoconversionalAnalysis.GAS_CONSTANT_KJ_PER_MOL_K * temperature),
        ) *
        reaction_function(model, conversion; m=m, n=n, p=p)
    for index in 1:(length(temperature_K) - 1)
        step = temperature_K[index + 1] - temperature_K[index]
        temperature = temperature_K[index]
        conversion = alpha[index]
        k1 = kinetic_rate(temperature, conversion) / beta
        k2 =
            kinetic_rate(
                temperature + step / 2,
                clamp(conversion + step * k1 / 2, 1.0e-8, 1 - 1.0e-8),
            ) / beta
        k3 =
            kinetic_rate(
                temperature + step / 2,
                clamp(conversion + step * k2 / 2, 1.0e-8, 1 - 1.0e-8),
            ) / beta
        k4 =
            kinetic_rate(
                temperature + step, clamp(conversion + step * k3, 1.0e-8, 1 - 1.0e-8)
            ) / beta
        alpha[index + 1] = clamp(
            conversion + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6, 1.0e-8, 1 - 1.0e-8
        )
    end
    exact_rate = [
        kinetic_rate(temperature, conversion) for
        (temperature, conversion) in zip(temperature_K, alpha)
    ]
    observed_rate = if noisy
        [
            rate * exp(
                0.015 * (
                    2 * mod(1_103_515_245 * index + 12_345, 2_147_483_648) / 2_147_483_648 - 1
                ),
            ) for (index, rate) in enumerate(exact_rate)
        ]
    else
        exact_rate
    end
    if inject_nonpositive
        observed_rate[findfirst(>(0.2), alpha)] = -1.0e-6
    end
    dalpha_dT = observed_rate ./ beta
    time_min = (temperature_K .- first(temperature_K)) ./ beta
    reconstructed = cumulative_trapezoid(temperature_K, dalpha_dT)
    residual = reconstructed .- (alpha .- first(alpha))
    row_count = length(alpha)
    preprocessing = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    ).preprocessing
    diagnostics = PreprocessingDiagnostics(
        row_count,
        row_count,
        row_count,
        row_count,
        0,
        0,
        0,
        2,
        10,
        10,
        100.0,
        50.0,
        beta,
        minimum(alpha),
        maximum(alpha),
        0,
        0,
        count(<(0), observed_rate),
        sqrt(mean(abs2, residual)),
        maximum(abs, residual),
        String[],
    )
    return ProcessedExperiment(
        identifier,
        "synthetic",
        identifier,
        repeat("d", 64),
        composition,
        collect(1:row_count),
        temperature_K,
        copy(time_min),
        time_min,
        100 .- 50 .* alpha,
        1 .- 0.5 .* alpha,
        alpha,
        copy(alpha),
        dalpha_dT,
        observed_rate,
        fill(beta, row_count),
        reconstructed,
        residual,
        first(temperature_K),
        last(temperature_K),
        preprocessing,
        "synthetic_m6_config",
        diagnostics,
    )
end

function kinetic_parameter(result::KineticTripletResult, name::Symbol)
    return only(
        parameter.estimate for parameter in result.parameters if parameter.name == name
    )
end

@testset "reaction-function registry and domains" begin
    registry = reaction_model_registry()
    @test length(registry) == 19
    @test allunique(getproperty.(registry, :name))
    @test reaction_function(:f1, 0.5) == 0.5
    @test reaction_function(:p2, 0.25) == 1.0
    @test reaction_function(:r2, 0.75) == 1.0
    @test reaction_function(:d4, 0.5) > 0
    @test reaction_function(:sestak_berggren_2, 0.4; m=0.3, n=1.2) > 0
    @test reaction_function(:sestak_berggren_3, 0.4; m=0.3, n=1.2, p=0.2) > 0
    @test_throws DomainError reaction_function(:f1, 0.0)
    @test_throws DomainError reaction_function(:f1, 1.0)
    @test_throws ArgumentError reaction_function(:unknown, 0.5)
end

@testset "exact multi-rate kinetic-triplet recovery" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    experiments = [
        synthetic_kinetic_experiment(;
            identifier="synthetic_kinetic_rate$(lpad(round(Int, rate), 2, '0'))",
            heating_rate_K_per_min=rate,
        ) for rate in (5.0, 10.0, 15.0, 20.0)
    ]
    result = fit_kinetic_triplet(experiments, :sestak_berggren_2, config)
    reversed_result = fit_kinetic_triplet(reverse(experiments), :sestak_berggren_2, config)

    @test abs(kinetic_parameter(result, :activation_energy_kJ_per_mol) - 145.0) < 1.0e-5
    @test abs(kinetic_parameter(result, :log_preexponential_per_min) - log(1.0e9)) < 1.0e-5
    @test abs(kinetic_parameter(result, :m) - 0.25) < 1.0e-5
    @test abs(kinetic_parameter(result, :n) - 1.10) < 1.0e-5
    @test result.diagnostics.log_rate_rmse < 1.0e-5
    @test result.diagnostics.cross_validation_log_rate_rmse < 1.0e-5
    @test result.experiment_ids == sort(result.experiment_ids)
    @test result.analysis_fingerprint == reversed_result.analysis_fingerprint
    @test kinetic_parameter(result, :activation_energy_kJ_per_mol) ≈
        kinetic_parameter(reversed_result, :activation_energy_kJ_per_mol) atol = 1.0e-10
    point = 300
    predicted = predict_rate(result, result.temperature_K[point], result.alpha[point])
    @test predicted ≈ result.observed_rate_per_min[point] rtol = 2.0e-6
    lower, central, upper = predict_rate_confidence_interval(
        result, result.temperature_K[point], result.alpha[point]
    )
    @test lower <= central <= upper
end

@testset "noisy model selection, visible exclusions, and failures" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    experiments = [
        synthetic_kinetic_experiment(;
            identifier="noisy_kinetic_rate$(lpad(round(Int, rate), 2, '0'))",
            heating_rate_K_per_min=rate,
            noisy=true,
            inject_nonpositive=rate == 5.0,
        ) for rate in (5.0, 10.0, 15.0, 20.0)
    ]
    comparison = compare_reaction_models(experiments, config)
    reversed_comparison = compare_reaction_models(reverse(experiments), config)
    selected = only(
        result for
        result in comparison.results if result.model.name == comparison.selected_model
    )

    @test comparison.selected_model == :sestak_berggren_2
    @test comparison.criterion_minimum_model == :sestak_berggren_2
    @test selected.diagnostics.excluded_nonpositive_rate_count == 1
    @test selected.diagnostics.cross_validation_log_rate_rmse < 0.05
    @test abs(kinetic_parameter(selected, :activation_energy_kJ_per_mol) - 145.0) < 2.0
    @test comparison.analysis_fingerprint == reversed_comparison.analysis_fingerprint
    @test comparison.criterion_values ≈ reversed_comparison.criterion_values atol = 1.0e-8

    mixed = synthetic_kinetic_experiment(;
        identifier="mixed_kinetic", composition=Composition(0.75)
    )
    @test_throws KineticModelError fit_kinetic_triplet(
        [experiments[1], experiments[2], mixed], :f1, config
    )
    @test_throws KineticModelError fit_kinetic_triplet(experiments[1:2], :f1, config)
    @test_throws KineticModelError fit_kinetic_triplet(experiments, :unknown, config)
end

@testset "compensation remains a visible single-rate diagnostic" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    experiment = synthetic_kinetic_experiment(;
        identifier="synthetic_f1", model=:f1, m=0.0, n=1.0, noisy=true
    )
    compensation = analyze_compensation(experiment, 145.0, config)
    first_order_index = findfirst(==(:f1), compensation.model_names)
    @test !isnothing(first_order_index)
    @test abs(
        compensation.model_activation_energies_kJ_per_mol[first_order_index] - 145.0
    ) < 2.0
    @test compensation.regression.r_squared > 0.95
    @test compensation.estimated_preexponential_per_min > 0
    @test compensation.confidence_lower_preexponential_per_min <=
        compensation.estimated_preexponential_per_min <=
        compensation.confidence_upper_preexponential_per_min
    @test any(contains("not evidence"), compensation.warnings)
end
