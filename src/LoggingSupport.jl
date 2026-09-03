"""
    make_logger(config; io=stderr) -> ConsoleLogger

Create the project logger at the configured minimum level. Julia logging key-value pairs
are preserved as structured metadata by the logging backend.
"""
function make_logger(config::AnalysisConfig; io::IO=stderr)
    return ConsoleLogger(io, config.project.log_level)
end

"""
    with_project_logger(f, config; io=stderr)

Run `f` with the logger selected by `config`.
"""
function with_project_logger(f::Function, config::AnalysisConfig; io::IO=stderr)
    return with_logger(make_logger(config; io=io)) do
        return f()
    end
end
