function synthetic_preprocessing_experiment(; segment_available=true, perturb_mass=false)
    temperature_celsius = collect(100.0:2.0:750.0)
    temperature_K = temperature_celsius .+ 273.15
    time_min = (temperature_celsius .- first(temperature_celsius)) ./ 10.0
    true_alpha = (temperature_celsius .- 150.0) ./ 550.0
    mass_percent = 100.0 .- 50.0 .* true_alpha
    perturb_mass && (mass_percent[findfirst(==(300.0), temperature_celsius)] += 1.0)
    segments = segment_available ? ifelse.(temperature_celsius .< 120.0, 1, 2) : nothing

    return Experiment(;
        id="synthetic_dynamic",
        role=:calibration,
        kind=:dynamic,
        source_file="synthetic.mat",
        source_variable="synthetic",
        source_sha256=repeat("a", 64),
        sample_id="synthetic",
        composition=Composition(0.5),
        temperature_K=vcat(temperature_K, NaN),
        time_min=vcat(time_min, NaN),
        mass_percent=vcat(mass_percent, NaN),
        segment=isnothing(segments) ? nothing : vcat(segments, missing),
        nominal_heating_rate_K_per_min=10.0,
        atmosphere="N2",
    )
end

function preprocessing_variant(
    original::PreprocessingConfig;
    profile=original.profile,
    invalid_row_policy=original.invalid_row_policy,
    segment_policy=original.segment_policy,
    rebase_time=original.rebase_time,
    reference_mass_method=original.reference_mass_method,
    reference_half_window_K=original.reference_half_window_K,
    smoothing_method=original.smoothing_method,
    smoothing_half_window_K=original.smoothing_half_window_K,
    local_polynomial_degree=original.local_polynomial_degree,
    derivative_method=original.derivative_method,
    heating_rate_method=original.heating_rate_method,
    monotonic_conversion_policy=original.monotonic_conversion_policy,
    analysis_conversion_range=original.analysis_conversion_range,
    reconstruction_rmse_tolerance=original.reconstruction_rmse_tolerance,
)
    return PreprocessingConfig(
        profile,
        invalid_row_policy,
        segment_policy,
        rebase_time,
        reference_mass_method,
        reference_half_window_K,
        smoothing_method,
        smoothing_half_window_K,
        local_polynomial_degree,
        derivative_method,
        heating_rate_method,
        monotonic_conversion_policy,
        analysis_conversion_range,
        reconstruction_rmse_tolerance,
    )
end

@testset "preprocessing numerical utilities" begin
    x = [0.0, 0.2, 0.9, 1.7, 3.0, 4.4, 6.0]
    y = @. 2.0 + 3.0 * x - 0.5 * x^2
    @test finite_difference_derivative(x, y) ≈ @.(3.0 - x) atol = 1.0e-12

    cubic = @. 1.0 - 2.0 * x + 0.4 * x^2 + 0.1 * x^3
    fitted, derivative = local_polynomial_estimate(x, cubic, x; half_window=10.0, degree=3)
    @test fitted ≈ cubic atol = 1.0e-11
    @test derivative ≈ @.(-2.0 + 0.8 * x + 0.3 * x^2) atol = 1.0e-10

    @test cumulative_trapezoid([0.0, 1.0, 3.0], [0.0, 2.0, 6.0]) == [0.0, 1.0, 9.0]
    @test_throws DimensionMismatch finite_difference_derivative([1.0, 2.0], [1.0])
    @test_throws ArgumentError finite_difference_derivative(
        [1.0, 1.0, 2.0], [1.0, 2.0, 3.0]
    )
    @test_throws ArgumentError local_polynomial_estimate(x, cubic, x; half_window=-1.0)
end

@testset "synthetic preprocessing pipeline" begin
    root = pkgdir(IsoconversionalAnalysis)
    config = load_config(joinpath(root, "config", "analysis_defaults.toml"))
    experiment = synthetic_preprocessing_experiment()
    original_temperature = copy(experiment.temperature_K)
    original_mass = copy(experiment.mass_percent)

    processed = preprocess(experiment, config)
    diagnostics = processed.diagnostics
    expected_alpha = (processed.temperature_K .- (150.0 + 273.15)) ./ 550.0

    @test processed.source_id == experiment.id
    @test first(processed.source_row_indices) == 26
    @test last(processed.source_row_indices) == 301
    @test first(processed.temperature_K) == 150.0 + 273.15
    @test last(processed.temperature_K) == 700.0 + 273.15
    @test first(processed.time_min) == 0.0
    @test first(processed.raw_time_min) == 5.0
    @test processed.alpha ≈ expected_alpha atol = 1.0e-12
    @test processed.analysis_alpha ≈ expected_alpha atol = 1.0e-10
    @test processed.dalpha_dT_K_inv ≈ fill(1 / 550, length(expected_alpha)) atol = 1.0e-10
    @test processed.dalpha_dt_min_inv ≈ fill(10 / 550, length(expected_alpha)) atol =
        1.0e-10
    @test processed.heating_rate_K_per_min ≈ fill(10.0, length(expected_alpha)) atol =
        1.0e-10
    @test processed.reconstructed_alpha ≈ expected_alpha atol = 1.0e-10
    @test maximum(abs, processed.reconstruction_residual) < 1.0e-10
    @test diagnostics.input_row_count == length(experiment.temperature_K)
    @test diagnostics.dropped_invalid_row_count == 1
    @test diagnostics.selected_segment == 2
    @test diagnostics.alpha_outside_unit_interval_count == 0
    @test diagnostics.conversion_reversal_count == 0
    @test diagnostics.negative_derivative_count == 0
    @test any(contains("dropped 1 non-finite"), diagnostics.warnings)
    @test length(processed.config_fingerprint) == 64
    @test processed.config_fingerprint ==
        preprocessing_fingerprint(config.preprocessing, 150.0 + 273.15, 700.0 + 273.15)
    @test processed.config_fingerprint !=
        preprocessing_fingerprint(config.preprocessing, 120.0 + 273.15, 700.0 + 273.15)
    @test isequal(experiment.temperature_K, original_temperature)
    @test isequal(experiment.mass_percent, original_mass)

    missing_segment = synthetic_preprocessing_experiment(; segment_available=false)
    fallback = preprocess(missing_segment, config)
    @test isnothing(fallback.diagnostics.selected_segment)
    @test any(contains("recorded segment is unavailable"), fallback.diagnostics.warnings)

    strict_invalid = preprocessing_variant(config.preprocessing; invalid_row_policy=:error)
    strict_config = AnalysisConfig(
        config.source_path,
        config.project,
        config.units,
        config.conversion,
        strict_invalid,
        config.isoconversional,
        config.deconvolution,
        config.reaction_models,
    )
    @test_throws PreprocessingError preprocess(experiment, strict_config)
    @test_throws PreprocessingError preprocess(
        experiment, config; final_temperature_celsius=800.0
    )
end

@testset "preprocessing candidate methods and visible diagnostics" begin
    root = pkgdir(IsoconversionalAnalysis)
    config = load_config(joinpath(root, "config", "analysis_defaults.toml"))
    experiment = synthetic_preprocessing_experiment(; perturb_mass=true)

    raw_finite_difference = preprocessing_variant(
        config.preprocessing;
        profile="raw_fd",
        reference_mass_method=:linear_interpolation,
        smoothing_method=:none,
        derivative_method=:finite_difference,
        heating_rate_method=:global_linear,
    )
    candidate_config = AnalysisConfig(
        config.source_path,
        config.project,
        config.units,
        config.conversion,
        raw_finite_difference,
        config.isoconversional,
        config.deconvolution,
        config.reaction_models,
    )
    candidate = preprocess(experiment, candidate_config)
    @test candidate.alpha == candidate.analysis_alpha
    @test candidate.diagnostics.conversion_reversal_count > 0
    @test candidate.diagnostics.negative_derivative_count > 0
    @test any(contains("decrease"), candidate.diagnostics.warnings)
    @test any(contains("derivative"), candidate.diagnostics.warnings)

    monotonic_error = preprocessing_variant(
        raw_finite_difference; monotonic_conversion_policy=:error
    )
    error_config = AnalysisConfig(
        config.source_path,
        config.project,
        config.units,
        config.conversion,
        monotonic_error,
        config.isoconversional,
        config.deconvolution,
        config.reaction_models,
    )
    @test_throws PreprocessingError preprocess(experiment, error_config)
end
