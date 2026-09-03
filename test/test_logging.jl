@testset "logging" begin
    config_path = joinpath(
        pkgdir(IsoconversionalAnalysis), "config", "analysis_defaults.toml"
    )
    config = load_config(config_path)
    buffer = IOBuffer()

    result = with_project_logger(config; io=buffer) do
        @debug "debug_event" run_id = "hidden"
        @info "configuration_loaded" source = config.source_path
        42
    end

    output = String(take!(buffer))
    @test result == 42
    @test !occursin("debug_event", output)
    @test occursin("configuration_loaded", output)
    @test occursin(config.source_path, output)
end
