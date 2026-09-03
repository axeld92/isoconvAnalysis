@testset "configuration" begin
    project_root = pkgdir(IsoconversionalAnalysis)
    config_path = joinpath(project_root, "config", "analysis_defaults.toml")
    config = load_config(config_path)

    @test config.source_path == abspath(config_path)
    @test config.project.output_directory == joinpath(project_root, "results")
    @test config.project.log_level == Logging.Info
    @test config.units.measured_temperature == "°C"
    @test config.units.thermodynamic_temperature == "K"
    @test config.units.activation_energy == "kJ/mol"
    @test config.conversion.initial_temperature_celsius == 150.0
    @test config.conversion.final_temperature_celsius == 700.0
    @test config.conversion.initial_sensitivity_celsius == [120.0]
    @test config.conversion.final_sensitivity_celsius == [650.0, 750.0]
    @test config.conversion.include_highest_common_temperature
    @test config.preprocessing.profile == "m3_default_v1"
    @test config.preprocessing.invalid_row_policy == :drop_nonfinite
    @test config.preprocessing.segment_policy == :recorded_ramp
    @test config.preprocessing.reference_mass_method == :local_robust_linear
    @test config.preprocessing.reference_half_window_K == 2.0
    @test config.preprocessing.smoothing_method == :local_polynomial
    @test config.preprocessing.smoothing_half_window_K == 5.0
    @test config.preprocessing.local_polynomial_degree == 3
    @test config.preprocessing.derivative_method == :local_polynomial
    @test config.preprocessing.heating_rate_method == :local_polynomial
    @test config.preprocessing.monotonic_conversion_policy == :diagnose
    @test config.preprocessing.analysis_conversion_range == (0.05, 0.95)
    @test config.preprocessing.reconstruction_rmse_tolerance == 0.01
    @test length(config.isoconversional.alpha_grid) == 17
    @test config.isoconversional.alpha_grid[1:2] == [0.10, 0.15]
    @test config.isoconversional.interpolation_policy == :first_upward_crossing
    @test config.isoconversional.confidence_level == 0.95
    @test config.isoconversional.minimum_experiments == 3
    @test config.isoconversional.friedman_minimum_rate_per_min == 1.0e-10
    @test config.isoconversional.minimum_r_squared_warning == 0.90
    @test config.isoconversional.advanced_vyazovkin.delta_alpha == 0.02
    @test config.isoconversional.advanced_vyazovkin.integration_points == 41
    @test config.isoconversional.advanced_vyazovkin.energy_bounds_kJ_per_mol ==
        (20.0, 500.0)
    @test config.deconvolution.peak_counts == [2, 3, 4]
    @test config.deconvolution.selection_criterion == :bic
    @test config.deconvolution.fit_conversion_range == (0.05, 0.95)
    @test config.deconvolution.maximum_fit_points == 600
    @test config.deconvolution.minimum_center_separation_K == 8.0
    @test config.deconvolution.skew_bounds == (-1.5, 1.5)
    @test config.deconvolution.width_bounds_K == (5.0, 250.0)
    @test config.deconvolution.maximum_height_factor == 3.0
    @test config.deconvolution.minimum_component_area_fraction == 0.02
    @test config.deconvolution.multistart_count == 5
    @test config.deconvolution.joint_multistart_count == 2
    @test config.deconvolution.optimizer_max_iterations == 600
    @test config.deconvolution.optimizer_gradient_tolerance == 1.0e-7
    @test config.deconvolution.optimizer_function_tolerance == 1.0e-12
    @test config.deconvolution.identifiability_condition_warning == 1.0e10
    @test config.deconvolution.boundary_fraction_tolerance == 0.01
    @test config.deconvolution.confidence_level == 0.95
    @test length(config.reaction_models.models) == 19
    @test first(config.reaction_models.models) == :f1
    @test last(config.reaction_models.models) == :sestak_berggren_3
    @test length(config.reaction_models.compensation_models) == 17
    @test config.reaction_models.fit_conversion_range == (0.10, 0.90)
    @test config.reaction_models.maximum_points_per_experiment == 250
    @test config.reaction_models.minimum_rate_per_min == 1.0e-10
    @test config.reaction_models.selection_criterion == :bic
    @test config.reaction_models.energy_bounds_kJ_per_mol == (20.0, 500.0)
    @test config.reaction_models.sestak_berggren_exponent_bounds == (-5.0, 5.0)
    @test config.reaction_models.minimum_experiments == 3
    @test config.reaction_models.confidence_level == 0.95
    @test config.reaction_models.maximum_cross_validation_log_rmse == 0.25

    mktempdir() do unrelated_directory
        from_another_directory = cd(unrelated_directory) do
            load_config(config_path)
        end
        @test from_another_directory.project.output_directory ==
            joinpath(project_root, "results")
    end

    missing_path = joinpath(project_root, "config", "missing.toml")
    @test_throws ConfigurationError load_config(missing_path)

    mktempdir() do directory
        invalid_path = joinpath(directory, "invalid.toml")
        default_text = read(config_path, String)
        write(
            invalid_path,
            replace(default_text, "log_level = \"info\"" => "log_level = \"verbose\""),
        )
        error = try
            load_config(invalid_path)
            nothing
        catch caught
            caught
        end
        @test error isa ConfigurationError
        @test occursin("project.log_level", sprint(showerror, error))

        invalid_window_path = joinpath(directory, "invalid_window.toml")
        write(
            invalid_window_path,
            replace(
                default_text,
                "smoothing_half_window_kelvin = 5.0" => "smoothing_half_window_kelvin = -5.0",
            ),
        )
        @test_throws ConfigurationError load_config(invalid_window_path)

        invalid_method_path = joinpath(directory, "invalid_method.toml")
        write(
            invalid_method_path,
            replace(
                default_text,
                "derivative_method = \"local_polynomial\"" => "derivative_method = \"mystery\"",
            ),
        )
        @test_throws ConfigurationError load_config(invalid_method_path)

        invalid_alpha_path = joinpath(directory, "invalid_alpha.toml")
        write(
            invalid_alpha_path,
            replace(default_text, "alpha_grid = [0.10, 0.15" => "alpha_grid = [0.01, 0.15"),
        )
        @test_throws ConfigurationError load_config(invalid_alpha_path)

        invalid_r_squared_path = joinpath(directory, "invalid_r_squared.toml")
        write(
            invalid_r_squared_path,
            replace(
                default_text,
                "minimum_r_squared_warning = 0.90" => "minimum_r_squared_warning = 1.1",
            ),
        )
        @test_throws ConfigurationError load_config(invalid_r_squared_path)

        invalid_bounds_path = joinpath(directory, "invalid_bounds.toml")
        write(
            invalid_bounds_path,
            replace(
                default_text,
                "energy_bounds_kilojoule_per_mole = [20.0, 500.0]" => "energy_bounds_kilojoule_per_mole = [500.0, 20.0]",
            ),
        )
        @test_throws ConfigurationError load_config(invalid_bounds_path)

        invalid_peak_counts_path = joinpath(directory, "invalid_peak_counts.toml")
        write(
            invalid_peak_counts_path,
            replace(default_text, "peak_counts = [2, 3, 4]" => "peak_counts = [2, 4, 3]"),
        )
        @test_throws ConfigurationError load_config(invalid_peak_counts_path)

        invalid_skew_bounds_path = joinpath(directory, "invalid_skew_bounds.toml")
        write(
            invalid_skew_bounds_path,
            replace(
                default_text, "skew_bounds = [-1.5, 1.5]" => "skew_bounds = [-2.0, 1.5]"
            ),
        )
        @test_throws ConfigurationError load_config(invalid_skew_bounds_path)

        invalid_area_fraction_path = joinpath(directory, "invalid_area_fraction.toml")
        write(
            invalid_area_fraction_path,
            replace(
                default_text,
                "minimum_component_area_fraction = 0.02" => "minimum_component_area_fraction = 0.6",
            ),
        )
        @test_throws ConfigurationError load_config(invalid_area_fraction_path)

        invalid_reaction_model_path = joinpath(directory, "invalid_reaction_model.toml")
        write(
            invalid_reaction_model_path,
            replace(default_text, "\"sestak_berggren_3\"]" => "\"mystery_model\"]"),
        )
        @test_throws ConfigurationError load_config(invalid_reaction_model_path)

        invalid_reaction_range_path = joinpath(directory, "invalid_reaction_range.toml")
        write(
            invalid_reaction_range_path,
            replace(
                default_text,
                "fit_conversion_range = [0.10, 0.90]" => "fit_conversion_range = [0.90, 0.10]",
            ),
        )
        @test_throws ConfigurationError load_config(invalid_reaction_range_path)

        malformed_path = joinpath(directory, "malformed.toml")
        write(malformed_path, "[project\nlog_level = \"info\"")
        @test_throws ConfigurationError load_config(malformed_path)
    end
end
