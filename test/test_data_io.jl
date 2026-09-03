using MAT
using TOML

include(joinpath(@__DIR__, "fixtures", "representative_data.jl"))

function fixture_catalog(directory::AbstractString; checksum=nothing)
    mat_path = write_representative_mat(joinpath(directory, "representative.mat"))
    actual_checksum = source_sha256(mat_path)
    template = read(joinpath(@__DIR__, "fixtures", "representative_catalog.toml"), String)
    catalog_text = replace(template, "__SHA256__" => something(checksum, actual_checksum))
    catalog_path = joinpath(directory, "catalog.toml")
    write(catalog_path, catalog_text)
    return catalog_path, mat_path
end

@testset "explicit MAT and catalog ingestion" begin
    mktempdir() do directory
        catalog_path, mat_path = fixture_catalog(directory)

        @test list_mat_variables(mat_path) ==
            ["dynamic_fixture", "not_a_matrix", "ramp_hold_fixture"]
        @test isequal(
            load_mat_variable(mat_path, "dynamic_fixture"), REPRESENTATIVE_DYNAMIC_MATRIX
        )
        @test_throws DataImportError load_mat_variable(mat_path, "missing")
        @test_throws DataImportError load_mat_variable(mat_path, "not_a_matrix")

        catalog = load_dataset_catalog(catalog_path)
        @test catalog.schema_version == 1
        @test length(catalog.sources) == 2
        @test sum(length(source.experiments) for source in catalog.sources) == 2

        experiments = load_experiments(catalog)
        @test getfield.(experiments, :id) ==
            ["fixture_dynamic_run", "fixture_ramp_hold_run"]
        dynamic = experiments[1]
        @test length(dynamic.temperature_K) == 4
        @test dynamic.temperature_K[1] == 298.15
        @test isnan(dynamic.temperature_K[end])
        @test dynamic.segment === nothing
        @test dynamic.composition.waste_tire == 0.25
        @test trailing_invalid_row_count(dynamic) == 1
        @test any(contains("missing_optional_segment"), dynamic.import_warnings)

        ramp_hold = load_experiment(catalog, "fixture_ramp_hold_run")
        @test ramp_hold.role == :validation
        @test ramp_hold.segment[1:4] == [1, 3, 3, 4]
        @test ismissing(ramp_hold.segment[end])
        @test ramp_hold.nominal_hold_temperature_celsius == 120.0
        @test getfield.(load_experiments(catalog; roles=[:validation]), :id) ==
            [ramp_hold.id]
        @test_throws DataImportError load_experiment(catalog, "unknown")

        audits = audit_catalog(catalog)
        @test length(audits) == 2
        @test audits[1].invalid_core_row_count == 1
        @test audits[2].measured_heating_rate_K_per_min ≈ 20.0

        markdown_path = joinpath(directory, "audit.md")
        toml_path = joinpath(directory, "audit.toml")
        outputs = write_audit_reports(
            catalog, audits; markdown_path=markdown_path, toml_path=toml_path
        )
        @test outputs.markdown == markdown_path
        @test occursin("Selected experiments: 2", read(markdown_path, String))
        machine_audit = TOML.parsefile(toml_path)
        @test machine_audit["summary"]["experiment_count"] == 2
        @test length(machine_audit["experiments"]) == 2

        first_markdown = read(markdown_path, String)
        first_toml = read(toml_path, String)
        write_audit_reports(
            catalog, audits; markdown_path=markdown_path, toml_path=toml_path
        )
        @test read(markdown_path, String) == first_markdown
        @test read(toml_path, String) == first_toml
    end

    mktempdir() do directory
        bad_catalog_path, _ = fixture_catalog(directory; checksum=repeat("0", 64))
        bad_catalog = load_dataset_catalog(bad_catalog_path)
        @test_throws DataImportError load_experiments(bad_catalog)
    end
end
