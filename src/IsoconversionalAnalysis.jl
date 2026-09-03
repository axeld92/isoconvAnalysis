"""
    IsoconversionalAnalysis

Reproducible tools for thermogravimetric and isoconversional kinetic analysis. The package
exposes validated configuration, MAT ingestion, raw-data audit, preprocessing,
isoconversional-analysis, Fraser–Suzuki deconvolution, and kinetic-triplet APIs through M6.
"""
module IsoconversionalAnalysis

using Dates
using Distributions
using LinearAlgebra
using Logging
using MAT
using Optim
using Printf
using SHA
using Statistics
using TOML

include("Types.jl")
include("Configuration.jl")
include("LoggingSupport.jl")
include("DataIO.jl")
include("DataAudit.jl")
include("Preprocessing.jl")
include("Isoconversional.jl")
include("PeakModels.jl")
include("Deconvolution.jl")
include("ReactionModels.jl")

export AdvancedVyazovkinConfig,
    AdvancedVyazovkinDiagnostics,
    AnalysisConfig,
    Composition,
    CompensationResult,
    ConfigurationError,
    ConversionConfig,
    ConversionSample,
    DataImportError,
    DatasetCatalog,
    DatasetConfigError,
    DatasetSourceSpec,
    DeconvolutionConfig,
    DeconvolutionDiagnostics,
    DeconvolutionError,
    DeconvolutionResult,
    Experiment,
    ExperimentAudit,
    ExperimentSpec,
    ExperimentValidationError,
    FraserSuzukiPeak,
    HeatingRateValidation,
    IsoconversionalConfig,
    IsoconversionalError,
    IsoconversionalPointDiagnostics,
    IsoconversionalResult,
    JointDeconvolutionDiagnostics,
    JointDeconvolutionResult,
    KineticModelComparison,
    KineticModelError,
    KineticParameterEstimate,
    KineticTripletDiagnostics,
    KineticTripletResult,
    LinearRegressionDiagnostics,
    PeakCountComparison,
    PeakParameterUncertainty,
    PreprocessingConfig,
    PreprocessingDiagnostics,
    PreprocessingError,
    ProcessedExperiment,
    ProjectConfig,
    ReactionModelConfig,
    ReactionModelSpec,
    UnitConfig,
    ValidationIssue,
    analyze,
    analyze_compensation,
    audit_catalog,
    audit_experiment,
    compare_peak_counts,
    compare_reaction_models,
    cumulative_trapezoid,
    finite_difference_derivative,
    fit_deconvolution,
    fit_joint_deconvolution,
    fit_kinetic_triplet,
    fraser_suzuki,
    fraser_suzuki_mixture,
    interpolate_at_conversion,
    list_mat_variables,
    load_config,
    load_dataset_catalog,
    load_experiment,
    load_experiments,
    load_mat_variable,
    local_polynomial_estimate,
    make_logger,
    pairwise_vyazovkin_objective,
    pairwise_vyazovkin_variance,
    predict_rate,
    predict_rate_confidence_interval,
    preprocess,
    preprocessing_fingerprint,
    reaction_function,
    reaction_model_registry,
    source_sha256,
    trailing_invalid_row_count,
    valid_row_mask,
    validate_experiment,
    with_project_logger,
    write_audit_reports

end
