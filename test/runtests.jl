using IsoconversionalAnalysis
using Logging
using Statistics
using Test

include(joinpath(@__DIR__, "..", "scripts", "audit_legacy_comparison.jl"))

@testset "IsoconversionalAnalysis" begin
    include("test_configuration.jl")
    include("test_logging.jl")
    include("test_data_io.jl")
    include("test_thermal_observations.jl")
    include("test_thermal_delimited.jl")
    include("test_thermal_tg_conversion.jl")
    include("test_thermal_conversion_support.jl")
    include("test_thermal_friedman_support.jl")
    include("test_thermal_vyazovkin_support.jl")
    include("test_thermal_ramp_support.jl")
    include("test_thermal_activation_energy.jl")
    include("test_experiment.jl")
    include("test_preprocessing.jl")
    include("test_isoconversional.jl")
    include("test_deconvolution.jl")
    include("test_reaction_models.jl")
    include("test_legacy_integration.jl")
    include("test_legacy_comparison.jl")
    include("test_quality.jl")
end
