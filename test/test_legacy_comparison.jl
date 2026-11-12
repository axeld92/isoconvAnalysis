using Dates
using SHA
using TOML

@testset "M9 corrected legacy comparison" begin
    audit = LegacyComparisonAudit
    project_root = audit.PROJECT_ROOT

    experiments = audit.collect_experiment_evidence(project_root)
    artifacts = audit.collect_artifact_evidence(project_root)
    @test length(experiments) == 20
    @test all(row -> row.alias_exact, experiments)
    @test all(row -> row.corrected_exact_contiguous_subset, experiments)
    @test all(row -> row.corrected_columns == 4, experiments)
    @test all(
        row -> row.corrected_end_row - row.corrected_start_row + 1 == row.corrected_rows,
        experiments,
    )
    @test minimum(row.corrected_start_temperature_celsius for row in experiments) > 109
    @test maximum(row.corrected_start_temperature_celsius for row in experiments) < 111
    @test minimum(row.corrected_end_temperature_celsius for row in experiments) > 790
    @test maximum(row.corrected_end_temperature_celsius for row in experiments) < 800

    @test length(audit.METHOD_COMPARISONS) == 8
    @test only(filter(row -> row.area == "Raw data ingestion", audit.METHOD_COMPARISONS)).status ==
        "exact"
    @test length(audit.INTENTIONAL_DEVIATIONS) == 13
    @test length(unique(row.id for row in audit.INTENTIONAL_DEVIATIONS)) == 13
    @test length(artifacts) == 18
    @test count(artifact -> artifact.tier == "saved_result_workspace", artifacts) == 10
    @test count(artifact -> artifact.tier == "publication_workspace", artifacts) == 2
    @test all(artifact -> occursin(r"^[0-9a-f]{64}$", artifact.sha256), artifacts)

    report = audit.render_report(experiments, artifacts; generated_on=Date(2026, 9, 5))
    @test startswith(report, "# M9 corrected Julia-versus-legacy comparison")
    @test occursin("20/20 canonical variables are exactly equal", report)
    @test occursin("All 20/20 arrays are exact contiguous copies", report)
    @test occursin("`LD04` | Friedman regression", report)
    @test occursin(
        "legacy comparison and intentional-deviation record are complete", report
    )
    @test !occursin(project_root, report)

    mktempdir() do directory
        report_path = joinpath(directory, "legacy_comparison.md")
        write(report_path, report)
        machine = audit.render_machine_evidence(
            project_root, experiments, artifacts, report_path; generated_on=Date(2026, 9, 5)
        )
        parsed = TOML.parse(machine)
        @test parsed["schema_version"] == 1
        @test parsed["report_sha256"] == audit.sha256_file(report_path)
        @test parsed["numerical_parity_scope"] == "raw_dynamic_import_only"
        @test !parsed["legacy_result_workspaces_are_numerical_oracles"]
        @test parsed["exact_raw_alias_count"] == 20
        @test parsed["exact_corrected_contiguous_subset_count"] == 20
        @test length(parsed["experiment"]) == 20
        @test length(parsed["method_comparison"]) == 8
        @test length(parsed["intentional_deviation"]) == 13
        @test length(parsed["artifact"]) == 18
        @test occursin(r"^[0-9a-f]{64}$", parsed["analysis_fingerprint"])
    end

    invalid = copy(experiments)
    invalid[1] = merge(invalid[1], (alias_exact=false,))
    @test_throws ArgumentError audit.validate_evidence(invalid, artifacts)
end
