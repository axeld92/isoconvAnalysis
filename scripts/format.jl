using JuliaFormatter

function main(args::Vector{String})
    project_root = normpath(joinpath(@__DIR__, ".."))
    targets = [
        joinpath(project_root, "src"),
        joinpath(project_root, "test"),
        joinpath(project_root, "scripts"),
        joinpath(project_root, "docs", "make.jl"),
    ]
    check_only = "--check" in args
    formatted = true

    for target in targets
        formatted &= JuliaFormatter.format(target; overwrite=(!check_only), verbose=true)
    end

    if check_only && !formatted
        error("Julia formatting check failed; run scripts/format.jl to apply formatting")
    end
end

main(ARGS)
