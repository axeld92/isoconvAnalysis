using IsoconversionalAnalysis
using Logging
using Statistics
using Test

@testset "IsoconversionalAnalysis" begin
    include("test_configuration.jl")
    include("test_logging.jl")
    include("test_data_io.jl")
    include("test_experiment.jl")
    include("test_preprocessing.jl")
    include("test_isoconversional.jl")
    include("test_deconvolution.jl")
    include("test_reaction_models.jl")
    include("test_legacy_integration.jl")
    include("test_quality.jl")
end
