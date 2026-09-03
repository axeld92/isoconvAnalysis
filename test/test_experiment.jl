@testset "Experiment validation" begin
    @test Composition(0.25) == Composition(0.25, 0.75)
    @test_throws ArgumentError Composition(0.6, 0.6)
    @test_throws ArgumentError Composition(-0.1)

    common = (
        id="validation_fixture",
        role=:calibration,
        kind=:dynamic,
        source_file="fixture.mat",
        source_variable="fixture",
        source_sha256=repeat("0", 64),
        sample_id="fixture",
        composition=Composition(0.5),
        temperature_K=[300.0, 310.0, NaN],
        time_min=[0.0, 1.0, NaN],
        mass_percent=[100.0, 90.0, NaN],
        nominal_heating_rate_K_per_min=10.0,
        atmosphere="N2",
    )
    experiment = Experiment(; common...)
    issues = validate_experiment(experiment)
    @test trailing_invalid_row_count(experiment) == 1
    @test count(valid_row_mask(experiment)) == 2
    @test any(issue -> issue.code == :trailing_invalid_rows, issues)
    @test any(contains("trailing_invalid_rows"), experiment.import_warnings)

    @test_throws ExperimentValidationError Experiment(; merge(common, (time_min=[0.0],))...)
    @test_throws ExperimentValidationError Experiment(;
        merge(common, (nominal_heating_rate_K_per_min=-1.0,))...
    )
    @test_throws ExperimentValidationError Experiment(;
        merge(common, (kind=:ramp_hold,))...
    )

    nonmonotonic = Experiment(; merge(common, (time_min=[0.0, 0.0, NaN],))...)
    @test any(contains("nonincreasing_time"), nonmonotonic.import_warnings)
end
