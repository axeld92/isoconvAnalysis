include(joinpath(@__DIR__, "..", "scripts", "audit_ictac_tg.jl"))

@testset "ICTAC TG assessment" begin
    @test ICTACTGAudit.program_checks([5, 5, 10])[1]["status"] == "not met"
    @test ICTACTGAudit.program_checks([5, 10, 20])[1]["status"] == "met"
    @test ICTACTGAudit.program_checks([5, 10, 20])[2]["status"] == "advisory"
    @test ICTACTGAudit.program_checks([2, 5, 20])[2]["status"] == "met"

    point(alpha, status, value) =
        Dict("alpha" => alpha, "status" => status, "energy_kJ_per_mol" => value)
    profile(case, points) = Dict("case" => case, "method" => "friedman", "points" => points)
    a = profile("reference", [point(0.1, "estimated", 100), point(0.2, "estimated", 200)])
    b = profile("sensitivity", [point(0.1, "estimated", 101), point(0.2, "warning", 500)])
    shift = ICTACTGAudit.profile_shift(a, b, [0.1, 0.9])
    @test shift["common_alpha"] == [0.1]
    @test shift["max_abs_shift_kJ_per_mol"] == 1

    config = ICTACTGAudit.load_toml(ICTACTGAudit.ASSESSMENT_PROTOCOL)
    tires = ICTACTGAudit.tire_cohort(config)
    @test tires.rates == [5, 10, 15, 20]
    @test all(c -> c.segment_label == "2" && c.mass_unit == :percent, tires.curves)
    @test length(unique(c.source_sha256 for c in tires.curves)) == 4
    @test all(
        m -> m["sample_mass_mg"] > 10 && !isempty(m["calibration_file"]), tires.metadata
    )
end
