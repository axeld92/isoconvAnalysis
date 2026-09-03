"""
Error raised when an analysis configuration is missing or invalid.
"""
struct ConfigurationError <: Exception
    path::String
    message::String
end

function Base.showerror(io::IO, error::ConfigurationError)
    return print(io, "invalid configuration at $(error.path): $(error.message)")
end

"""
Project-level output and logging settings.
"""
struct ProjectConfig
    output_directory::String
    log_level::LogLevel
end

"""
Canonical unit labels used at file and reporting boundaries.
"""
struct UnitConfig
    measured_temperature::String
    thermodynamic_temperature::String
    time::String
    mass::String
    heating_rate::String
    activation_energy::String
    conversion::String
    conversion_rate::String
end

"""
Primary conversion interval and the endpoint sensitivity settings.
"""
struct ConversionConfig
    initial_temperature_celsius::Float64
    final_temperature_celsius::Float64
    initial_sensitivity_celsius::Vector{Float64}
    final_sensitivity_celsius::Vector{Float64}
    include_highest_common_temperature::Bool
end

"""
Explicit numerical choices used to transform a raw experiment.

Temperature windows are half-widths in kelvin. `monotonic_conversion_policy = :diagnose`
records reversals without altering the data; no supported policy silently clips or repairs
conversion.
"""
struct PreprocessingConfig
    profile::String
    invalid_row_policy::Symbol
    segment_policy::Symbol
    rebase_time::Bool
    reference_mass_method::Symbol
    reference_half_window_K::Float64
    smoothing_method::Symbol
    smoothing_half_window_K::Float64
    local_polynomial_degree::Int
    derivative_method::Symbol
    heating_rate_method::Symbol
    monotonic_conversion_policy::Symbol
    analysis_conversion_range::Tuple{Float64,Float64}
    reconstruction_rmse_tolerance::Float64
end

"""
Numerical integration, bounded optimization, and diagnostic settings for the advanced
Vyazovkin method.
"""
struct AdvancedVyazovkinConfig
    delta_alpha::Float64
    integration_points::Int
    energy_bounds_kJ_per_mol::Tuple{Float64,Float64}
    optimizer_grid_points::Int
    optimizer_tolerance_kJ_per_mol::Float64
    optimizer_max_iterations::Int
    boundary_tolerance_kJ_per_mol::Float64
    flat_relative_tolerance::Float64
end

"""
Shared conversion grid, interpolation, uncertainty, and method-specific settings for M4
isoconversional calculations.
"""
struct IsoconversionalConfig
    alpha_grid::Vector{Float64}
    interpolation_policy::Symbol
    confidence_level::Float64
    minimum_experiments::Int
    friedman_minimum_rate_per_min::Float64
    minimum_r_squared_warning::Float64
    advanced_vyazovkin::AdvancedVyazovkinConfig
end

"""
Peak-count comparison, physical bounds, deterministic multistart, optimizer, and
identifiability settings for M5 Fraser–Suzuki deconvolution.
"""
struct DeconvolutionConfig
    peak_counts::Vector{Int}
    selection_criterion::Symbol
    fit_conversion_range::Tuple{Float64,Float64}
    maximum_fit_points::Int
    minimum_center_separation_K::Float64
    skew_bounds::Tuple{Float64,Float64}
    width_bounds_K::Tuple{Float64,Float64}
    maximum_height_factor::Float64
    minimum_component_area_fraction::Float64
    multistart_count::Int
    joint_multistart_count::Int
    optimizer_max_iterations::Int
    optimizer_gradient_tolerance::Float64
    optimizer_function_tolerance::Float64
    identifiability_condition_warning::Float64
    boundary_fraction_tolerance::Float64
    confidence_level::Float64
end

"""
Candidate reaction models, shared fit interval, physical plausibility bounds, uncertainty,
model-selection, and leave-one-heating-rate-out validation settings for M6.
"""
struct ReactionModelConfig
    models::Vector{Symbol}
    compensation_models::Vector{Symbol}
    fit_conversion_range::Tuple{Float64,Float64}
    maximum_points_per_experiment::Int
    minimum_rate_per_min::Float64
    selection_criterion::Symbol
    energy_bounds_kJ_per_mol::Tuple{Float64,Float64}
    sestak_berggren_exponent_bounds::Tuple{Float64,Float64}
    minimum_experiments::Int
    confidence_level::Float64
    identifiability_condition_warning::Float64
    parameter_correlation_warning::Float64
    ambiguity_criterion_delta::Float64
    maximum_cross_validation_log_rmse::Float64
end

"""
Validated configuration for an analysis run.
"""
struct AnalysisConfig
    source_path::String
    project::ProjectConfig
    units::UnitConfig
    conversion::ConversionConfig
    preprocessing::PreprocessingConfig
    isoconversional::IsoconversionalConfig
    deconvolution::DeconvolutionConfig
    reaction_models::ReactionModelConfig
end

const _FIXED_REACTION_MODELS = (
    :f1,
    :f2,
    :f3,
    :f4,
    :a2,
    :a3,
    :a4,
    :r2,
    :r3,
    :d1,
    :d2,
    :d3,
    :d4,
    :p2,
    :p3,
    :p4,
    :random_scission,
)
const _EMPIRICAL_REACTION_MODELS = (:sestak_berggren_2, :sestak_berggren_3)
const _SUPPORTED_REACTION_MODELS = (
    _FIXED_REACTION_MODELS..., _EMPIRICAL_REACTION_MODELS...
)

const _EXPECTED_UNITS = (
    measured_temperature="°C",
    thermodynamic_temperature="K",
    time="min",
    mass="mg",
    heating_rate="K/min",
    activation_energy="kJ/mol",
    conversion="1",
    conversion_rate="min^-1",
)

function _config_error(path::AbstractString, message::AbstractString)
    return throw(ConfigurationError(String(path), String(message)))
end

function _required_table(data::AbstractDict, key::String, path::String)
    haskey(data, key) || _config_error(path, "missing [$key] table")
    value = data[key]
    value isa AbstractDict || _config_error(path, "[$key] must be a TOML table")
    return value
end

function _required_string(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    value = table[key]
    value isa AbstractString || _config_error(path, "$section.$key must be a string")
    isempty(strip(value)) && _config_error(path, "$section.$key cannot be empty")
    return String(value)
end

function _required_number(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    value = table[key]
    value isa Real && !(value isa Bool) ||
        _config_error(path, "$section.$key must be a number")
    isfinite(value) || _config_error(path, "$section.$key must be finite")
    return Float64(value)
end

function _required_numbers(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    values = table[key]
    values isa AbstractVector || _config_error(path, "$section.$key must be an array")
    result = Float64[]
    for value in values
        value isa Real && !(value isa Bool) ||
            _config_error(path, "$section.$key must contain only numbers")
        isfinite(value) ||
            _config_error(path, "$section.$key must contain only finite values")
        push!(result, Float64(value))
    end
    isempty(result) && _config_error(path, "$section.$key cannot be empty")
    return result
end

function _required_strings(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    values = table[key]
    values isa AbstractVector || _config_error(path, "$section.$key must be an array")
    result = String[]
    for value in values
        value isa AbstractString ||
            _config_error(path, "$section.$key must contain only strings")
        isempty(strip(value)) && _config_error(path, "$section.$key cannot contain blanks")
        push!(result, String(value))
    end
    isempty(result) && _config_error(path, "$section.$key cannot be empty")
    return result
end

function _required_bool(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    value = table[key]
    value isa Bool || _config_error(path, "$section.$key must be true or false")
    return value
end

function _required_integer(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    value = table[key]
    value isa Integer && !(value isa Bool) ||
        _config_error(path, "$section.$key must be an integer")
    return Int(value)
end

function _required_integers(table::AbstractDict, key::String, section::String, path::String)
    haskey(table, key) || _config_error(path, "missing $section.$key")
    values = table[key]
    values isa AbstractVector || _config_error(path, "$section.$key must be an array")
    result = Int[]
    for value in values
        value isa Integer && !(value isa Bool) ||
            _config_error(path, "$section.$key must contain only integers")
        push!(result, Int(value))
    end
    isempty(result) && _config_error(path, "$section.$key cannot be empty")
    return result
end

function _required_symbol(
    table::AbstractDict, key::String, section::String, path::String, allowed::Tuple
)
    value = Symbol(_required_string(table, key, section, path))
    value in allowed || _config_error(
        path,
        "$section.$key must be one of $(join(string.(allowed), ", ")); received $value",
    )
    return value
end

function _parse_log_level(value::AbstractString, path::String)
    levels = Dict(
        "debug" => Logging.Debug,
        "info" => Logging.Info,
        "warn" => Logging.Warn,
        "error" => Logging.Error,
    )
    normalized = lowercase(strip(value))
    haskey(levels, normalized) ||
        _config_error(path, "project.log_level must be debug, info, warn, or error")
    return levels[normalized]
end

function _parse_units(table::AbstractDict, path::String)
    values = map(keys(_EXPECTED_UNITS)) do key
        value = _required_string(table, String(key), "units", path)
        expected = getproperty(_EXPECTED_UNITS, key)
        value == expected ||
            _config_error(path, "units.$key must be '$expected'; received '$value'")
        return value
    end
    return UnitConfig(values...)
end

function _parse_conversion(table::AbstractDict, path::String)
    initial = _required_number(table, "initial_temperature_celsius", "conversion", path)
    final = _required_number(table, "final_temperature_celsius", "conversion", path)
    initial < final || _config_error(
        path, "conversion initial temperature must be below final temperature"
    )

    initial_sensitivity = _required_numbers(
        table, "initial_sensitivity_celsius", "conversion", path
    )
    final_sensitivity = _required_numbers(
        table, "final_sensitivity_celsius", "conversion", path
    )
    all(initial_sensitivity .< final) || _config_error(
        path,
        "all conversion initial-temperature sensitivities must be below the primary final temperature",
    )
    all(final_sensitivity .> initial) || _config_error(
        path,
        "all conversion final-temperature sensitivities must exceed the primary initial temperature",
    )

    include_highest_common = _required_bool(
        table, "include_highest_common_temperature", "conversion", path
    )
    return ConversionConfig(
        initial, final, initial_sensitivity, final_sensitivity, include_highest_common
    )
end

function _parse_preprocessing(table::AbstractDict, path::String)
    section = "preprocessing"
    profile = _required_string(table, "profile", section, path)
    invalid_policy = _required_symbol(
        table, "invalid_row_policy", section, path, (:drop_nonfinite, :error)
    )
    segment_policy = _required_symbol(
        table, "segment_policy", section, path, (:recorded_ramp, :all_finite)
    )
    rebase_time = _required_bool(table, "rebase_time", section, path)
    reference_method = _required_symbol(
        table,
        "reference_mass_method",
        section,
        path,
        (:local_robust_linear, :linear_interpolation),
    )
    reference_window = _required_number(
        table, "reference_half_window_kelvin", section, path
    )
    reference_window > 0 ||
        _config_error(path, "$section.reference_half_window_kelvin must be positive")
    smoothing_method = _required_symbol(
        table, "smoothing_method", section, path, (:none, :local_polynomial)
    )
    smoothing_window = _required_number(
        table, "smoothing_half_window_kelvin", section, path
    )
    smoothing_window > 0 ||
        _config_error(path, "$section.smoothing_half_window_kelvin must be positive")
    degree = _required_integer(table, "local_polynomial_degree", section, path)
    1 <= degree <= 5 ||
        _config_error(path, "$section.local_polynomial_degree must be between 1 and 5")
    derivative_method = _required_symbol(
        table, "derivative_method", section, path, (:finite_difference, :local_polynomial)
    )
    heating_rate_method = _required_symbol(
        table, "heating_rate_method", section, path, (:global_linear, :local_polynomial)
    )
    monotonic_policy = _required_symbol(
        table, "monotonic_conversion_policy", section, path, (:diagnose, :error)
    )
    conversion_range = _required_numbers(table, "analysis_conversion_range", section, path)
    length(conversion_range) == 2 ||
        _config_error(path, "$section.analysis_conversion_range must have two values")
    0 <= conversion_range[1] < conversion_range[2] <= 1 || _config_error(
        path, "$section.analysis_conversion_range must satisfy 0 <= lower < upper <= 1"
    )
    reconstruction_tolerance = _required_number(
        table, "reconstruction_rmse_tolerance", section, path
    )
    reconstruction_tolerance > 0 ||
        _config_error(path, "$section.reconstruction_rmse_tolerance must be positive")

    return PreprocessingConfig(
        profile,
        invalid_policy,
        segment_policy,
        rebase_time,
        reference_method,
        reference_window,
        smoothing_method,
        smoothing_window,
        degree,
        derivative_method,
        heating_rate_method,
        monotonic_policy,
        (conversion_range[1], conversion_range[2]),
        reconstruction_tolerance,
    )
end

function _parse_advanced_vyazovkin(table::AbstractDict, path::String)
    section = "isoconversional.advanced_vyazovkin"
    delta_alpha = _required_number(table, "delta_alpha", section, path)
    0 < delta_alpha < 1 ||
        _config_error(path, "$section.delta_alpha must lie strictly between zero and one")
    integration_points = _required_integer(table, "integration_points", section, path)
    integration_points >= 3 ||
        _config_error(path, "$section.integration_points must be at least three")
    bounds = _required_numbers(table, "energy_bounds_kilojoule_per_mole", section, path)
    length(bounds) == 2 || _config_error(
        path, "$section.energy_bounds_kilojoule_per_mole must have two values"
    )
    0 < bounds[1] < bounds[2] || _config_error(
        path,
        "$section.energy_bounds_kilojoule_per_mole must contain increasing positive values",
    )
    grid_points = _required_integer(table, "optimizer_grid_points", section, path)
    grid_points >= 5 ||
        _config_error(path, "$section.optimizer_grid_points must be at least five")
    tolerance = _required_number(
        table, "optimizer_tolerance_kilojoule_per_mole", section, path
    )
    tolerance > 0 || _config_error(
        path, "$section.optimizer_tolerance_kilojoule_per_mole must be positive"
    )
    maximum_iterations = _required_integer(table, "optimizer_max_iterations", section, path)
    maximum_iterations > 0 ||
        _config_error(path, "$section.optimizer_max_iterations must be positive")
    boundary_tolerance = _required_number(
        table, "boundary_tolerance_kilojoule_per_mole", section, path
    )
    0 <= boundary_tolerance < (bounds[2] - bounds[1]) / 2 || _config_error(
        path,
        "$section.boundary_tolerance_kilojoule_per_mole must be nonnegative and smaller than half the energy range",
    )
    flat_tolerance = _required_number(table, "flat_relative_tolerance", section, path)
    flat_tolerance > 0 ||
        _config_error(path, "$section.flat_relative_tolerance must be positive")
    return AdvancedVyazovkinConfig(
        delta_alpha,
        integration_points,
        (bounds[1], bounds[2]),
        grid_points,
        tolerance,
        maximum_iterations,
        boundary_tolerance,
        flat_tolerance,
    )
end

function _parse_isoconversional(table::AbstractDict, path::String)
    section = "isoconversional"
    alpha_grid = _required_numbers(table, "alpha_grid", section, path)
    all(value -> 0 < value < 1, alpha_grid) || _config_error(
        path, "$section.alpha_grid values must lie strictly between zero and one"
    )
    all(>(0), diff(alpha_grid)) ||
        _config_error(path, "$section.alpha_grid must be strictly increasing")
    interpolation_policy = _required_symbol(
        table, "interpolation_policy", section, path, (:first_upward_crossing,)
    )
    confidence_level = _required_number(table, "confidence_level", section, path)
    0 < confidence_level < 1 || _config_error(
        path, "$section.confidence_level must lie strictly between zero and one"
    )
    minimum_experiments = _required_integer(table, "minimum_experiments", section, path)
    minimum_experiments >= 3 ||
        _config_error(path, "$section.minimum_experiments must be at least three")
    minimum_rate = _required_number(
        table, "friedman_minimum_rate_per_minute", section, path
    )
    minimum_rate > 0 ||
        _config_error(path, "$section.friedman_minimum_rate_per_minute must be positive")
    minimum_r_squared_warning = _required_number(
        table, "minimum_r_squared_warning", section, path
    )
    0 <= minimum_r_squared_warning <= 1 || _config_error(
        path, "$section.minimum_r_squared_warning must lie between zero and one"
    )
    advanced_table = _required_table(table, "advanced_vyazovkin", path)
    advanced = _parse_advanced_vyazovkin(advanced_table, path)
    minimum(alpha_grid) > advanced.delta_alpha || _config_error(
        path, "$section.alpha_grid must start above advanced_vyazovkin.delta_alpha"
    )
    return IsoconversionalConfig(
        alpha_grid,
        interpolation_policy,
        confidence_level,
        minimum_experiments,
        minimum_rate,
        minimum_r_squared_warning,
        advanced,
    )
end

function _parse_deconvolution(table::AbstractDict, path::String)
    section = "deconvolution"
    peak_counts = _required_integers(table, "peak_counts", section, path)
    issorted(peak_counts) &&
    allunique(peak_counts) &&
    all(count -> 1 <= count <= 6, peak_counts) || _config_error(
        path, "$section.peak_counts must be unique, increasing, and lie in 1:6"
    )
    selection_criterion = _required_symbol(
        table, "selection_criterion", section, path, (:bic, :aicc)
    )
    conversion_range = _required_numbers(table, "fit_conversion_range", section, path)
    length(conversion_range) == 2 ||
        _config_error(path, "$section.fit_conversion_range must have two values")
    0 <= conversion_range[1] < conversion_range[2] <= 1 || _config_error(
        path, "$section.fit_conversion_range must satisfy 0 <= lower < upper <= 1"
    )
    maximum_fit_points = _required_integer(table, "maximum_fit_points", section, path)
    maximum_fit_points >= 50 ||
        _config_error(path, "$section.maximum_fit_points must be at least 50")
    minimum_separation = _required_number(
        table, "minimum_center_separation_kelvin", section, path
    )
    minimum_separation > 0 ||
        _config_error(path, "$section.minimum_center_separation_kelvin must be positive")
    skew_bounds = _required_numbers(table, "skew_bounds", section, path)
    length(skew_bounds) == 2 && -2 < skew_bounds[1] < skew_bounds[2] < 2 || _config_error(
        path, "$section.skew_bounds must be increasing and lie within (-2, 2)"
    )
    width_bounds = _required_numbers(table, "width_bounds_kelvin", section, path)
    length(width_bounds) == 2 && 0 < width_bounds[1] < width_bounds[2] || _config_error(
        path, "$section.width_bounds_kelvin must contain increasing positive values"
    )
    maximum_height_factor = _required_number(table, "maximum_height_factor", section, path)
    maximum_height_factor > 1 ||
        _config_error(path, "$section.maximum_height_factor must exceed one")
    minimum_area_fraction = _required_number(
        table, "minimum_component_area_fraction", section, path
    )
    0 < minimum_area_fraction < 0.5 || _config_error(
        path, "$section.minimum_component_area_fraction must lie between zero and 0.5"
    )
    multistart_count = _required_integer(table, "multistart_count", section, path)
    multistart_count >= 1 ||
        _config_error(path, "$section.multistart_count must be positive")
    joint_multistart_count = _required_integer(
        table, "joint_multistart_count", section, path
    )
    joint_multistart_count >= 1 ||
        _config_error(path, "$section.joint_multistart_count must be positive")
    maximum_iterations = _required_integer(table, "optimizer_max_iterations", section, path)
    maximum_iterations >= 1 ||
        _config_error(path, "$section.optimizer_max_iterations must be positive")
    gradient_tolerance = _required_number(
        table, "optimizer_gradient_tolerance", section, path
    )
    gradient_tolerance > 0 ||
        _config_error(path, "$section.optimizer_gradient_tolerance must be positive")
    function_tolerance = _required_number(
        table, "optimizer_function_tolerance", section, path
    )
    function_tolerance > 0 ||
        _config_error(path, "$section.optimizer_function_tolerance must be positive")
    condition_warning = _required_number(
        table, "identifiability_condition_warning", section, path
    )
    condition_warning > 1 ||
        _config_error(path, "$section.identifiability_condition_warning must exceed one")
    boundary_fraction = _required_number(
        table, "boundary_fraction_tolerance", section, path
    )
    0 < boundary_fraction < 0.5 || _config_error(
        path, "$section.boundary_fraction_tolerance must lie between zero and 0.5"
    )
    confidence_level = _required_number(table, "confidence_level", section, path)
    0 < confidence_level < 1 ||
        _config_error(path, "$section.confidence_level must lie between zero and one")
    return DeconvolutionConfig(
        peak_counts,
        selection_criterion,
        (conversion_range[1], conversion_range[2]),
        maximum_fit_points,
        minimum_separation,
        (skew_bounds[1], skew_bounds[2]),
        (width_bounds[1], width_bounds[2]),
        maximum_height_factor,
        minimum_area_fraction,
        multistart_count,
        joint_multistart_count,
        maximum_iterations,
        gradient_tolerance,
        function_tolerance,
        condition_warning,
        boundary_fraction,
        confidence_level,
    )
end

function _parse_reaction_models(table::AbstractDict, path::String)
    section = "reaction_models"
    model_strings = _required_strings(table, "models", section, path)
    models = Symbol.(model_strings)
    allunique(models) || _config_error(path, "$section.models must be unique")
    all(model -> model in _SUPPORTED_REACTION_MODELS, models) || _config_error(
        path,
        "$section.models contains an unsupported model; allowed values are $(join(string.(_SUPPORTED_REACTION_MODELS), ", "))",
    )
    compensation_strings = _required_strings(table, "compensation_models", section, path)
    compensation_models = Symbol.(compensation_strings)
    allunique(compensation_models) ||
        _config_error(path, "$section.compensation_models must be unique")
    all(model -> model in _FIXED_REACTION_MODELS, compensation_models) || _config_error(
        path, "$section.compensation_models must contain only fixed registry models"
    )
    length(compensation_models) >= 3 || _config_error(
        path, "$section.compensation_models must contain at least three models"
    )
    conversion_range = _required_numbers(table, "fit_conversion_range", section, path)
    length(conversion_range) == 2 ||
        _config_error(path, "$section.fit_conversion_range must have two values")
    0 < conversion_range[1] < conversion_range[2] < 1 ||
        _config_error(path, "$section.fit_conversion_range must lie strictly inside (0, 1)")
    maximum_points = _required_integer(
        table, "maximum_points_per_experiment", section, path
    )
    maximum_points >= 20 ||
        _config_error(path, "$section.maximum_points_per_experiment must be at least 20")
    minimum_rate = _required_number(table, "minimum_rate_per_minute", section, path)
    minimum_rate > 0 ||
        _config_error(path, "$section.minimum_rate_per_minute must be positive")
    selection_criterion = _required_symbol(
        table, "selection_criterion", section, path, (:bic, :aicc)
    )
    energy_bounds = _required_numbers(
        table, "energy_bounds_kilojoule_per_mole", section, path
    )
    length(energy_bounds) == 2 && 0 < energy_bounds[1] < energy_bounds[2] || _config_error(
        path,
        "$section.energy_bounds_kilojoule_per_mole must contain two increasing positive values",
    )
    exponent_bounds = _required_numbers(
        table, "sestak_berggren_exponent_bounds", section, path
    )
    length(exponent_bounds) == 2 && exponent_bounds[1] < exponent_bounds[2] ||
        _config_error(path, "$section.sestak_berggren_exponent_bounds must be increasing")
    minimum_experiments = _required_integer(table, "minimum_experiments", section, path)
    minimum_experiments >= 3 ||
        _config_error(path, "$section.minimum_experiments must be at least three")
    confidence_level = _required_number(table, "confidence_level", section, path)
    0 < confidence_level < 1 || _config_error(
        path, "$section.confidence_level must lie strictly between zero and one"
    )
    condition_warning = _required_number(
        table, "identifiability_condition_warning", section, path
    )
    condition_warning > 1 ||
        _config_error(path, "$section.identifiability_condition_warning must exceed one")
    correlation_warning = _required_number(
        table, "parameter_correlation_warning", section, path
    )
    0 < correlation_warning < 1 || _config_error(
        path, "$section.parameter_correlation_warning must lie between zero and one"
    )
    ambiguity_delta = _required_number(table, "ambiguity_criterion_delta", section, path)
    ambiguity_delta >= 0 ||
        _config_error(path, "$section.ambiguity_criterion_delta must be nonnegative")
    maximum_cv_rmse = _required_number(
        table, "maximum_cross_validation_log_rmse", section, path
    )
    maximum_cv_rmse > 0 ||
        _config_error(path, "$section.maximum_cross_validation_log_rmse must be positive")
    return ReactionModelConfig(
        models,
        compensation_models,
        (conversion_range[1], conversion_range[2]),
        maximum_points,
        minimum_rate,
        selection_criterion,
        (energy_bounds[1], energy_bounds[2]),
        (exponent_bounds[1], exponent_bounds[2]),
        minimum_experiments,
        confidence_level,
        condition_warning,
        correlation_warning,
        ambiguity_delta,
        maximum_cv_rmse,
    )
end

"""
    load_config(path) -> AnalysisConfig

Parse and validate an analysis configuration from TOML. Relative output paths are
resolved relative to the configuration file, so behavior is independent of the current
working directory. Unit labels are currently restricted to the canonical boundary units
recorded in `docs/data_dictionary.md`; numerical code will convert temperatures to kelvin.
"""
function load_config(path::AbstractString)
    absolute_path = abspath(path)
    isfile(absolute_path) || _config_error(absolute_path, "file does not exist")

    data = try
        TOML.parsefile(absolute_path)
    catch error
        error isa TOML.ParserError || rethrow()
        _config_error(absolute_path, sprint(showerror, error))
    end

    project_table = _required_table(data, "project", absolute_path)
    units_table = _required_table(data, "units", absolute_path)
    conversion_table = _required_table(data, "conversion", absolute_path)
    preprocessing_table = _required_table(data, "preprocessing", absolute_path)
    isoconversional_table = _required_table(data, "isoconversional", absolute_path)
    deconvolution_table = _required_table(data, "deconvolution", absolute_path)
    reaction_model_table = _required_table(data, "reaction_models", absolute_path)

    raw_output = _required_string(
        project_table, "output_directory", "project", absolute_path
    )
    output_directory = if isabspath(raw_output)
        normpath(raw_output)
    else
        normpath(joinpath(dirname(absolute_path), raw_output))
    end
    log_level = _parse_log_level(
        _required_string(project_table, "log_level", "project", absolute_path),
        absolute_path,
    )

    return AnalysisConfig(
        absolute_path,
        ProjectConfig(output_directory, log_level),
        _parse_units(units_table, absolute_path),
        _parse_conversion(conversion_table, absolute_path),
        _parse_preprocessing(preprocessing_table, absolute_path),
        _parse_isoconversional(isoconversional_table, absolute_path),
        _parse_deconvolution(deconvolution_table, absolute_path),
        _parse_reaction_models(reaction_model_table, absolute_path),
    )
end
