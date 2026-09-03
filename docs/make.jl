using Documenter
using IsoconversionalAnalysis

makedocs(;
    sitename="IsoconversionalAnalysis.jl",
    format=Documenter.HTML(; edit_link=nothing, repolink=nothing),
    modules=[IsoconversionalAnalysis],
    source=joinpath(@__DIR__, "src"),
    build=joinpath(@__DIR__, "build"),
    clean=true,
    doctest=true,
    checkdocs=:exports,
    remotes=nothing,
    pages=[
        "Home" => "index.md",
        "Configuration and logging" => "configuration.md",
        "Data ingestion and audit" => "data_ingestion.md",
        "Preprocessing" => "preprocessing.md",
        "Isoconversional analysis" => "isoconversional.md",
        "Fraser–Suzuki deconvolution" => "deconvolution.md",
        "Kinetic triplets and reaction models" => "reaction_models.md",
        "Development" => "development.md",
    ],
)
