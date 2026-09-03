function synthetic_deconvolution_experiment(;
    identifier="synthetic_deconvolution",
    center_shift_K=0.0,
    noisy=false,
    composition=Composition(0.5),
)
    temperature_K = collect(430.0:1.0:970.0)
    generating_peaks = [
        FraserSuzukiPeak(0.0040, -0.20, 555.0 + center_shift_K, 75.0),
        FraserSuzukiPeak(0.0060, 0.16, 680.0 + center_shift_K, 95.0),
        FraserSuzukiPeak(0.0032, -0.10, 805.0 + center_shift_K, 115.0),
    ]
    raw_rate = fraser_suzuki_mixture(temperature_K, generating_peaks)
    normalization = last(cumulative_trapezoid(temperature_K, raw_rate))
    exact_rate = raw_rate ./ normalization
    peaks = [
        FraserSuzukiPeak(
            peak.height_K_inv / normalization, peak.skew, peak.center_K, peak.width_K
        ) for peak in generating_peaks
    ]
    observed_rate = if noisy
        [
            rate +
            2.0e-5 * (
                2 * mod(1_103_515_245 * index + 12_345, 2_147_483_648) / 2_147_483_648 - 1
            ) for (index, rate) in enumerate(exact_rate)
        ]
    else
        exact_rate
    end
    alpha = cumulative_trapezoid(temperature_K, exact_rate)
    beta = 10.0
    time_min = (temperature_K .- first(temperature_K)) ./ beta
    reconstructed = cumulative_trapezoid(temperature_K, observed_rate)
    residual = reconstructed .- alpha
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
        count(<(-1.0e-10), observed_rate),
        sqrt(sum(abs2, residual) / row_count),
        maximum(abs, residual),
        String[],
    )
    experiment = ProcessedExperiment(
        identifier,
        "synthetic",
        identifier,
        repeat("c", 64),
        composition,
        collect(1:row_count),
        temperature_K,
        copy(time_min),
        time_min,
        100 .- 50 .* alpha,
        1 .- 0.5 .* alpha,
        alpha,
        copy(alpha),
        observed_rate,
        beta .* observed_rate,
        fill(beta, row_count),
        reconstructed,
        residual,
        first(temperature_K),
        last(temperature_K),
        preprocessing,
        "synthetic_deconvolution_config",
        diagnostics,
    )
    return experiment, peaks
end

function deconvolution_configuration_variant(
    original::DeconvolutionConfig;
    peak_counts=original.peak_counts,
    selection_criterion=original.selection_criterion,
    fit_conversion_range=original.fit_conversion_range,
    maximum_fit_points=original.maximum_fit_points,
    minimum_center_separation_K=original.minimum_center_separation_K,
    skew_bounds=original.skew_bounds,
    width_bounds_K=original.width_bounds_K,
    maximum_height_factor=original.maximum_height_factor,
    minimum_component_area_fraction=original.minimum_component_area_fraction,
    multistart_count=original.multistart_count,
    joint_multistart_count=original.joint_multistart_count,
    optimizer_max_iterations=original.optimizer_max_iterations,
    optimizer_gradient_tolerance=original.optimizer_gradient_tolerance,
    optimizer_function_tolerance=original.optimizer_function_tolerance,
    identifiability_condition_warning=original.identifiability_condition_warning,
    boundary_fraction_tolerance=original.boundary_fraction_tolerance,
    confidence_level=original.confidence_level,
)
    return DeconvolutionConfig(
        peak_counts,
        selection_criterion,
        fit_conversion_range,
        maximum_fit_points,
        minimum_center_separation_K,
        skew_bounds,
        width_bounds_K,
        maximum_height_factor,
        minimum_component_area_fraction,
        multistart_count,
        joint_multistart_count,
        optimizer_max_iterations,
        optimizer_gradient_tolerance,
        optimizer_function_tolerance,
        identifiability_condition_warning,
        boundary_fraction_tolerance,
        confidence_level,
    )
end

@testset "Fraser–Suzuki kernel contracts" begin
    gaussian = FraserSuzukiPeak(2.0, 0.0, 600.0, 80.0)
    @test fraser_suzuki(600.0, gaussian) == 2.0
    @test fraser_suzuki(560.0, gaussian) ≈ 1.0 atol = 1.0e-14
    @test fraser_suzuki(640.0, gaussian) ≈ 1.0 atol = 1.0e-14
    almost_gaussian = FraserSuzukiPeak(2.0, 1.0e-9, 600.0, 80.0)
    @test fraser_suzuki(collect(500.0:10.0:700.0), almost_gaussian) ≈
        fraser_suzuki(collect(500.0:10.0:700.0), gaussian)

    skewed = FraserSuzukiPeak(1.0, 0.5, 600.0, 80.0)
    @test fraser_suzuki(500.0, skewed) == 0.0
    @test fraser_suzuki(600.0, skewed) == 1.0
    @test_throws ArgumentError FraserSuzukiPeak(-1.0, 0.0, 600.0, 80.0)
    @test_throws ArgumentError FraserSuzukiPeak(1.0, 0.0, 600.0, 0.0)
    @test_throws ArgumentError fraser_suzuki_mixture([500.0], FraserSuzukiPeak[])
end

@testset "known three-peak synthetic deconvolution" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    experiment, truth = synthetic_deconvolution_experiment()
    original_temperature = copy(experiment.temperature_K)
    original_rate = copy(experiment.dalpha_dT_K_inv)
    result = fit_deconvolution(experiment, 3, config)

    @test result.peak_count == 3
    @test issorted(getproperty.(result.peaks, :center_K))
    @test maximum(
        abs.(getproperty.(result.peaks, :center_K) .- getproperty.(truth, :center_K))
    ) < 4.0
    @test maximum(
        abs.(
            result.component_area_fractions .- begin
                areas = [
                    last(
                        cumulative_trapezoid(
                            result.temperature_K, fraser_suzuki(result.temperature_K, peak)
                        ),
                    ) for peak in truth
                ]
                areas ./ sum(areas)
            end,
        ),
    ) < 0.04
    @test result.diagnostics.root_mean_square_error < 5.0e-5
    @test result.diagnostics.relative_area_error < 0.01
    @test length(result.parameter_uncertainty) == 3
    @test size(result.parameter_covariance) == (12, 12)
    @test length(result.analysis_fingerprint) == 64
    @test experiment.temperature_K == original_temperature
    @test experiment.dalpha_dT_K_inv == original_rate
end

@testset "peak-count selection and visible input failures" begin
    config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    experiment, _ = synthetic_deconvolution_experiment(; noisy=true)
    comparison = compare_peak_counts(experiment, config)
    @test comparison.selected_peak_count == 3
    @test length(comparison.results) == 3
    @test minimum(comparison.criterion_deltas) == 0.0
    @test_throws DeconvolutionError fit_deconvolution(experiment, 5, config)
end

@testset "joint shared-skew recovery and experiment-order invariance" begin
    analysis_config = load_config(
        joinpath(pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml")
    )
    config = deconvolution_configuration_variant(
        analysis_config.deconvolution;
        multistart_count=3,
        joint_multistart_count=1,
        optimizer_max_iterations=400,
    )
    first_experiment, truth = synthetic_deconvolution_experiment(;
        identifier="synthetic_rate05", center_shift_K=-4.0, noisy=true
    )
    second_experiment, _ = synthetic_deconvolution_experiment(;
        identifier="synthetic_rate20", center_shift_K=4.0, noisy=true
    )
    experiments = [first_experiment, second_experiment]
    result = fit_joint_deconvolution(experiments, 3; configuration=config)
    reversed_result = fit_joint_deconvolution(reverse(experiments), 3; configuration=config)

    @test result.experiment_ids == ["synthetic_rate05", "synthetic_rate20"]
    @test maximum(abs.(result.shared_skews .- getproperty.(truth, :skew))) < 0.08
    @test all(
        experiment_result -> issorted(getproperty.(experiment_result.peaks, :center_K)),
        result.experiment_results,
    )
    @test result.analysis_fingerprint == reversed_result.analysis_fingerprint
    @test result.shared_skews ≈ reversed_result.shared_skews atol = 1.0e-8
    @test result.diagnostics.relative_objective_increase < 0.15

    mixed_experiment, _ = synthetic_deconvolution_experiment(;
        identifier="mixed", composition=Composition(0.75)
    )
    @test_throws DeconvolutionError fit_joint_deconvolution(
        [first_experiment, mixed_experiment], 3; configuration=config
    )
    @test_throws DeconvolutionError fit_joint_deconvolution(
        experiments, 5; configuration=config
    )
    @test_throws DeconvolutionError fit_joint_deconvolution(
        experiments, 3; configuration=config, independent_results=result.experiment_results
    )
end
