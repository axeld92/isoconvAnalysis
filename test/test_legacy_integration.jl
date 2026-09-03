@testset "selected legacy dataset integration" begin
    project_root = pkgdir(IsoconversionalAnalysis)
    catalog = load_dataset_catalog(joinpath(project_root, "config", "datasets.toml"))
    experiments = load_experiments(catalog)

    @test length(experiments) == 35
    @test count(experiment -> experiment.role == :calibration, experiments) == 20
    @test count(experiment -> experiment.role == :validation, experiments) == 15
    @test all(experiment -> trailing_invalid_row_count(experiment) == 1, experiments)
    @test all(experiment -> count(valid_row_mask(experiment)) >= 5999, experiments)
    @test all(
        experiment -> experiment.source_sha256 == source_sha256(experiment.source_file),
        experiments,
    )

    no_segment_ids = sort!([
        experiment.id for experiment in experiments if isnothing(experiment.segment)
    ])
    @test no_segment_ids ==
        ["dynamic_wt000_rate05", "dynamic_wt000_rate10", "dynamic_wt000_rate15"]
    @test all(
        experiment ->
            isnothing(experiment.raw_export_file) || isfile(experiment.raw_export_file),
        experiments,
    )

    audits = audit_experiment.(experiments)
    dynamic_audits = filter(audit -> audit.kind == :dynamic, audits)
    ramp_hold_audits = filter(audit -> audit.kind == :ramp_hold, audits)
    @test maximum(
        abs(audit.measured_heating_rate_K_per_min - audit.nominal_heating_rate_K_per_min)
        for audit in dynamic_audits
    ) < 0.5
    @test all(audit -> audit.time_nonincreasing_step_count == 0, audits)
    @test all(
        audit ->
            audit.median_purge_flow_mL_per_min == 80.0 &&
            audit.median_protective_flow_mL_per_min == 20.0,
        dynamic_audits,
    )
    @test all(
        audit ->
            audit.median_purge_flow_mL_per_min == 20.0 &&
            audit.median_protective_flow_mL_per_min == 20.0,
        ramp_hold_audits,
    )
    @test all(audit -> audit.segment_values == [1, 2, 3, 4], ramp_hold_audits)
    @test all(
        audit -> 20.6 < audit.measured_heating_rate_K_per_min < 21.0, ramp_hold_audits
    )

    source_variables = Dict(
        source.id => list_mat_variables(source.source_file) for source in catalog.sources
    )
    @test "ExpDatdrifttest2" in source_variables["dynamic_2021"]
    @test length(source_variables["ramp_hold_2021"]) == 15
end

@testset "M4–M6 methods run on one complete real calibration group" begin
    project_root = pkgdir(IsoconversionalAnalysis)
    config = load_config(joinpath(project_root, "config", "analysis_defaults.toml"))
    catalog = load_dataset_catalog(joinpath(project_root, "config", "datasets.toml"))
    experiments = load_experiments(catalog; roles=[:calibration])
    group = filter(experiments) do experiment
        return experiment.composition.waste_tire == 0.5
    end
    processed = preprocess.(group, Ref(config))

    @test length(processed) == 4
    for method in (:friedman, :kas, :fwo, :starink, :advanced_vyazovkin)
        result = analyze(method, processed, config)
        @test length(result.alpha) == 17
        @test all(!ismissing, result.activation_energy_kJ_per_mol)
        @test all(isfinite, Float64.(result.activation_energy_kJ_per_mol))
        @test result.experiment_ids == sort(result.experiment_ids)
        @test all(
            point.status in (:ok, :warning, :boundary_solution, :multimodal_objective) for
            point in result.point_diagnostics
        )
        if method == :friedman
            @test any(point -> point.status == :warning, result.point_diagnostics)
            @test any(
                point -> any(contains("R²"), point.warnings), result.point_diagnostics
            )
        end
    end

    representative = only(
        filter(result -> result.source_id == "dynamic_wt050_rate10", processed)
    )
    comparison = compare_peak_counts(representative, config)
    @test comparison.criterion_minimum_peak_count == 4
    @test comparison.selected_peak_count == 3
    @test comparison.status == :constraint_filtered
    @test comparison.structurally_eligible == [false, true, false]
    @test comparison.results[2].diagnostics.root_mean_square_error < 2.0e-4
    @test comparison.results[2].diagnostics.relative_area_error < 0.01

    kinetic_comparison = compare_reaction_models(processed, config)
    @test length(kinetic_comparison.results) == 19
    @test kinetic_comparison.criterion_minimum_model == :sestak_berggren_3
    @test kinetic_comparison.selected_model == :sestak_berggren_2
    @test kinetic_comparison.status == :predictive_warning
    selected_kinetic = only(
        result for result in kinetic_comparison.results if
        result.model.name == kinetic_comparison.selected_model
    )
    @test selected_kinetic.diagnostics.cross_validation_log_rate_rmse >
        config.reaction_models.maximum_cross_validation_log_rmse
    @test kinetic_comparison.analysis_fingerprint ==
        compare_reaction_models(reverse(processed), config).analysis_fingerprint
end

@testset "all dynamic runs preprocess with the accepted M3 profile" begin
    project_root = pkgdir(IsoconversionalAnalysis)
    config = load_config(joinpath(project_root, "config", "analysis_defaults.toml"))
    catalog = load_dataset_catalog(joinpath(project_root, "config", "datasets.toml"))
    experiments = load_experiments(catalog; roles=[:calibration])
    processed = preprocess.(experiments, Ref(config))

    @test length(processed) == 20
    @test all(result -> result.initial_temperature_K == 423.15, processed)
    @test all(result -> result.final_temperature_K == 973.15, processed)
    @test all(result -> all(isfinite, result.dalpha_dT_K_inv), processed)
    @test all(result -> all(>(0), result.heating_rate_K_per_min), processed)
    @test maximum(result.diagnostics.reconstruction_rmse for result in processed) < 0.001
    @test maximum(
        abs(
            result.diagnostics.measured_heating_rate_K_per_min -
            experiment.nominal_heating_rate_K_per_min,
        ) for (result, experiment) in zip(processed, experiments)
    ) < 0.5
    @test all(
        result -> result.source_row_indices == sort(result.source_row_indices), processed
    )
    @test count(
        result ->
            any(contains("recorded segment is unavailable"), result.diagnostics.warnings),
        processed,
    ) == 3
    @test any(result -> result.diagnostics.alpha_outside_unit_interval_count > 0, processed)
    @test length(unique(result.config_fingerprint for result in processed)) == 1
end
