function synthetic_kinetic_experiment(
    heating_rate_K_per_min::Real;
    activation_energy_kJ_per_mol=150.0,
    preexponential_per_min=1.0e12,
    noisy=false,
    composition=Composition(0.5),
    identifier=nothing,
)
    beta = Float64(heating_rate_K_per_min)
    temperature_K = collect(350.0:0.25:1000.0)
    time_min = (temperature_K .- first(temperature_K)) ./ beta
    arrhenius = @. exp(
        -activation_energy_kJ_per_mol /
        (IsoconversionalAnalysis.GAS_CONSTANT_KJ_PER_MOL_K * temperature_K),
    )
    thermal_integral = cumulative_trapezoid(temperature_K, arrhenius)
    alpha = @. 1 - exp(-(preexponential_per_min / beta) * thermal_integral)
    analysis_alpha = if noisy
        [
            value +
            value *
            (1 - value) *
            (2.0e-3 * sin(0.071 * index) + 1.0e-3 * cos(0.037 * index)) for
            (index, value) in enumerate(alpha)
        ]
    else
        copy(alpha)
    end
    exact_rate = @. preexponential_per_min * arrhenius * (1 - alpha)
    dalpha_dt = if noisy
        [rate * (1 + 0.01 * sin(0.053 * index)) for (index, rate) in enumerate(exact_rate)]
    else
        exact_rate
    end
    dalpha_dT = dalpha_dt ./ beta
    reconstructed = cumulative_trapezoid(temperature_K, dalpha_dT; initial=alpha[1])
    residual = reconstructed .- analysis_alpha
    row_count = length(alpha)
    configuration = load_config(
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
        count(<(-1.0e-8), diff(analysis_alpha)),
        count(<(-1.0e-10), dalpha_dT),
        sqrt(sum(abs2, residual) / row_count),
        maximum(abs, residual),
        String[],
    )
    rate_label = lpad(round(Int, beta), 2, '0')
    source_id = isnothing(identifier) ? "synthetic_rate$rate_label" : String(identifier)
    return ProcessedExperiment(
        source_id,
        "synthetic",
        "rate$rate_label",
        repeat("b", 64),
        composition,
        collect(1:row_count),
        temperature_K,
        copy(time_min),
        time_min,
        (@. 100 - 50 * alpha),
        (@. 1 - 0.5 * alpha),
        alpha,
        analysis_alpha,
        dalpha_dT,
        dalpha_dt,
        fill(beta, row_count),
        reconstructed,
        residual,
        first(temperature_K),
        last(temperature_K),
        configuration,
        "synthetic_first_order",
        diagnostics,
    )
end

function isoconversional_configuration_variant(
    original::IsoconversionalConfig;
    alpha_grid=original.alpha_grid,
    interpolation_policy=original.interpolation_policy,
    confidence_level=original.confidence_level,
    minimum_experiments=original.minimum_experiments,
    friedman_minimum_rate_per_min=original.friedman_minimum_rate_per_min,
    minimum_r_squared_warning=original.minimum_r_squared_warning,
    advanced_vyazovkin=original.advanced_vyazovkin,
)
    return IsoconversionalConfig(
        alpha_grid,
        interpolation_policy,
        confidence_level,
        minimum_experiments,
        friedman_minimum_rate_per_min,
        minimum_r_squared_warning,
        advanced_vyazovkin,
    )
end

@testset "shared isoconversional interpolation and pairwise objective" begin
    experiment = synthetic_kinetic_experiment(10.0)
    sample = interpolate_at_conversion(experiment, 0.5)
    @test sample.alpha == 0.5
    @test 350.0 < sample.temperature_K < 1000.0
    @test sample.time_min > 0
    @test sample.dalpha_dt_min_inv > 0
    @test sample.crossing_count == 1
    limited_conversion = synthetic_kinetic_experiment(10.0; preexponential_per_min=1.0e4)
    @test_throws IsoconversionalError interpolate_at_conversion(limited_conversion, 0.99)
    @test_throws IsoconversionalError interpolate_at_conversion(
        experiment, 0.5; policy=:spline
    )

    @test pairwise_vyazovkin_objective(zeros(3)) ≈ 6.0
    @test pairwise_vyazovkin_objective(zeros(4)) ≈ 12.0
    @test pairwise_vyazovkin_objective(zeros(5)) ≈ 20.0
    @test pairwise_vyazovkin_variance(zeros(5)) == 0.0
    @test pairwise_vyazovkin_objective([0.0]) == Inf
    @test pairwise_vyazovkin_variance([0.0]) == Inf
    @test IsoconversionalAnalysis._crossing_indices([0.0, 0.5, 0.4], 0.5) == [1]
    @test IsoconversionalAnalysis._crossing_indices([0.5, 0.6, 0.7], 0.5) == [1]
end

@testset "known-energy first-order synthetic recovery" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    alpha_grid = [0.2, 0.4, 0.6, 0.8]
    experiments = synthetic_kinetic_experiment.([5.0, 10.0, 15.0, 20.0])
    methods = (:friedman, :kas, :fwo, :starink, :advanced_vyazovkin)
    results = Dict(
        method =>
            analyze(method, experiments, alpha_grid; configuration=config.isoconversional)
        for method in methods
    )

    @test all(
        result -> all(!ismissing, result.activation_energy_kJ_per_mol), values(results)
    )
    @test maximum(
        abs.(Float64.(results[:friedman].activation_energy_kJ_per_mol) .- 150.0)
    ) < 0.5
    @test maximum(
        abs.(Float64.(results[:advanced_vyazovkin].activation_energy_kJ_per_mol) .- 150.0)
    ) < 0.75
    for method in (:kas, :fwo, :starink)
        @test maximum(
            abs.(Float64.(results[method].activation_energy_kJ_per_mol) .- 150.0)
        ) < 8.0
    end
    @test all(
        point -> point.regression.r_squared > 0.999, results[:friedman].point_diagnostics
    )
    @test all(
        point ->
            point.advanced_vyazovkin.objective_minimum + 1.0e-8 >=
            point.advanced_vyazovkin.objective_baseline,
        results[:advanced_vyazovkin].point_diagnostics,
    )
    @test results[:friedman].interval_method ==
        "two-sided Student-t interval for transformed OLS slope"
    @test results[:advanced_vyazovkin].interval_method ==
        "Vyazovkin–Wight pairwise-variance Fisher interval"
    @test all(result -> length(result.analysis_fingerprint) == 64, values(results))

    for method in methods
        reversed_result = analyze(
            method, reverse(experiments), alpha_grid; configuration=config.isoconversional
        )
        @test Float64.(reversed_result.activation_energy_kJ_per_mol) ≈
            Float64.(results[method].activation_energy_kJ_per_mol) atol = 1.0e-8
        @test reversed_result.analysis_fingerprint == results[method].analysis_fingerprint
    end
    advanced_result = results[:advanced_vyazovkin]
    @test all(
        !diagnostics.advanced_vyazovkin.fisher_lower_truncated &&
            !diagnostics.advanced_vyazovkin.fisher_upper_truncated for
        diagnostics in advanced_result.point_diagnostics
    )
    @test all(
        lower <= energy <= upper for (lower, energy, upper) in zip(
            advanced_result.confidence_lower_kJ_per_mol,
            advanced_result.activation_energy_kJ_per_mol,
            advanced_result.confidence_upper_kJ_per_mol,
        )
    )
end

@testset "noisy synthetic recovery and visible failures" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    alpha_grid = [0.2, 0.4, 0.6, 0.8]
    noisy = [
        synthetic_kinetic_experiment(rate; noisy=true) for rate in (5.0, 10.0, 15.0, 20.0)
    ]
    for method in (:friedman, :kas, :fwo, :starink, :advanced_vyazovkin)
        result = analyze(method, noisy, alpha_grid; configuration=config.isoconversional)
        estimates = Float64.(result.activation_energy_kJ_per_mol)
        @test maximum(abs.(estimates .- 150.0)) < 15.0
    end

    exact = synthetic_kinetic_experiment.([5.0, 10.0, 15.0, 20.0])
    @test_throws IsoconversionalError analyze(
        :mystery, exact, alpha_grid; configuration=config.isoconversional
    )
    @test_throws IsoconversionalError analyze(
        :kas, exact[1:2], alpha_grid; configuration=config.isoconversional
    )
    mixed_composition = copy(exact)
    mixed_composition[end] = synthetic_kinetic_experiment(
        20.0; composition=Composition(0.75)
    )
    @test_throws IsoconversionalError analyze(
        :kas, mixed_composition, alpha_grid; configuration=config.isoconversional
    )
    @test_throws IsoconversionalError analyze(
        :kas, exact, [0.4, 0.3]; configuration=config.isoconversional
    )

    indistinguishable = [
        synthetic_kinetic_experiment(10.0; identifier="same_rate_$index") for index in 1:4
    ]
    invalid = analyze(:kas, indistinguishable, [0.5]; configuration=config.isoconversional)
    @test ismissing(only(invalid.activation_energy_kJ_per_mol))
    @test only(invalid.point_diagnostics).status == :invalid
end

@testset "advanced Vyazovkin boundary diagnostics" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    experiments = [
        synthetic_kinetic_experiment(
            rate; activation_energy_kJ_per_mol=550.0, preexponential_per_min=1.0e30
        ) for rate in (5.0, 10.0, 15.0, 20.0)
    ]
    result = analyze(
        :advanced_vyazovkin, experiments, [0.5]; configuration=config.isoconversional
    )
    @test result.point_diagnostics[1].status == :boundary_solution
    @test result.point_diagnostics[1].advanced_vyazovkin.boundary_solution
    @test any(contains("boundary"), result.point_diagnostics[1].warnings)
end
