const GAS_CONSTANT_KJ_PER_MOL_K = 0.00831446261815324
const _SUPPORTED_ISOCONVERSIONAL_METHODS = (
    :friedman, :kas, :fwo, :starink, :advanced_vyazovkin
)

"""
Raised when an isoconversional calculation violates an input, interpolation, or numerical
contract.
"""
struct IsoconversionalError <: Exception
    method::Symbol
    message::String
end

function Base.showerror(io::IO, error::IsoconversionalError)
    return print(io, "$(error.method) analysis failed: $(error.message)")
end

"""
One acquisition-order interpolation at a requested conversion. `crossing_count` exposes
multiple upward crossings rather than silently treating a non-monotonic curve as invertible.
"""
struct ConversionSample
    alpha::Float64
    temperature_K::Float64
    time_min::Float64
    dalpha_dT_K_inv::Float64
    dalpha_dt_min_inv::Float64
    heating_rate_K_per_min::Float64
    crossing_count::Int
end

"""
Ordinary-least-squares diagnostics for one conversion value. The confidence interval stored
in the parent result is obtained by transforming the two-sided Student-t interval for
`slope`.
"""
struct LinearRegressionDiagnostics
    slope::Float64
    intercept::Float64
    slope_standard_error::Float64
    residual_standard_error::Float64
    r_squared::Float64
    degrees_of_freedom::Int
    x::Vector{Float64}
    y::Vector{Float64}
    fitted_y::Vector{Float64}
    residuals::Vector{Float64}
end

"""
Optimization and Fisher-interval diagnostics for one advanced Vyazovkin estimate. The energy
estimate minimizes the pairwise objective itself. `fisher_center_energy_kJ_per_mol` separately
records the minimum of the Vyazovkin–Wight pairwise variance used to construct the confidence
interval.
"""
struct AdvancedVyazovkinDiagnostics
    objective_minimum::Float64
    objective_baseline::Float64
    objective_evaluation_count::Int
    scan_local_minimum_count::Int
    relative_scan_range::Float64
    boundary_solution::Bool
    flat_objective::Bool
    fisher_variance_minimum::Float64
    fisher_center_energy_kJ_per_mol::Float64
    fisher_critical_value::Float64
    fisher_lower_truncated::Bool
    fisher_upper_truncated::Bool
end

"""
Per-conversion inputs, status, warnings, and method-specific diagnostics. A non-`:ok` status
is never silently dropped from the parent result.
"""
struct IsoconversionalPointDiagnostics
    alpha::Float64
    status::Symbol
    temperatures_K::Vector{Float64}
    global_heating_rates_K_per_min::Vector{Float64}
    conversion_rates_min_inv::Vector{Float64}
    crossing_counts::Vector{Int}
    regression::Union{Nothing,LinearRegressionDiagnostics}
    advanced_vyazovkin::Union{Nothing,AdvancedVyazovkinDiagnostics}
    warnings::Vector{String}
end

"""
Typed output of one isoconversional method over a common conversion grid. Missing estimates
are paired with a non-`:ok` point status. Interval construction and confidence level are
named explicitly, and provenance connects the result to every processed experiment and M4
configuration.
"""
struct IsoconversionalResult
    method::Symbol
    alpha::Vector{Float64}
    activation_energy_kJ_per_mol::Vector{Union{Missing,Float64}}
    confidence_lower_kJ_per_mol::Vector{Union{Missing,Float64}}
    confidence_upper_kJ_per_mol::Vector{Union{Missing,Float64}}
    confidence_level::Float64
    interval_method::String
    experiment_ids::Vector{String}
    preprocessing_fingerprints::Vector{String}
    configuration::IsoconversionalConfig
    point_diagnostics::Vector{IsoconversionalPointDiagnostics}
    analysis_fingerprint::String
end

function _isoconversional_error(method::Symbol, message::AbstractString)
    return throw(IsoconversionalError(method, String(message)))
end

function _crossing_indices(alpha::Vector{Float64}, target::Float64)
    indices = Int[]
    for index in 1:(length(alpha) - 1)
        lower = alpha[index]
        upper = alpha[index + 1]
        if upper > lower && (lower < target <= upper || (index == 1 && lower == target))
            push!(indices, index)
        end
    end
    return indices
end

"""
    interpolate_at_conversion(experiment, alpha; policy=:first_upward_crossing)

Interpolate temperature, time, conversion rates, and local heating rate at a common
conversion using acquisition-order upward crossings. Extrapolation is rejected. When noise
creates multiple crossings, the earliest is returned and the total is recorded in
`ConversionSample.crossing_count`.
"""
function interpolate_at_conversion(
    experiment::ProcessedExperiment, alpha::Real; policy::Symbol=:first_upward_crossing
)
    policy == :first_upward_crossing ||
        _isoconversional_error(:interpolation, "unsupported interpolation policy $policy")
    target = Float64(alpha)
    isfinite(target) && 0 <= target <= 1 || _isoconversional_error(
        :interpolation, "conversion must be finite and lie in [0, 1]"
    )
    crossings = _crossing_indices(experiment.analysis_alpha, target)
    isempty(crossings) && _isoconversional_error(
        :interpolation,
        "$(experiment.source_id) has no upward crossing at alpha=$target; extrapolation is disabled",
    )
    index = first(crossings)
    lower_alpha = experiment.analysis_alpha[index]
    upper_alpha = experiment.analysis_alpha[index + 1]
    fraction = (target - lower_alpha) / (upper_alpha - lower_alpha)
    interpolate(values) = values[index] + fraction * (values[index + 1] - values[index])
    return ConversionSample(
        target,
        interpolate(experiment.temperature_K),
        interpolate(experiment.time_min),
        interpolate(experiment.dalpha_dT_K_inv),
        interpolate(experiment.dalpha_dt_min_inv),
        interpolate(experiment.heating_rate_K_per_min),
        length(crossings),
    )
end

function _validate_analysis_inputs(
    method::Symbol,
    experiments::AbstractVector{<:ProcessedExperiment},
    alpha_grid::Vector{Float64},
    configuration::IsoconversionalConfig,
)
    method in _SUPPORTED_ISOCONVERSIONAL_METHODS || _isoconversional_error(
        method,
        "supported methods are $(join(string.(_SUPPORTED_ISOCONVERSIONAL_METHODS), ", "))",
    )
    length(experiments) >= configuration.minimum_experiments || _isoconversional_error(
        method,
        "at least $(configuration.minimum_experiments) experiments are required; received $(length(experiments))",
    )
    length(unique(experiment.source_id for experiment in experiments)) ==
    length(experiments) ||
        _isoconversional_error(method, "experiment identifiers must be unique")
    isempty(alpha_grid) && _isoconversional_error(method, "alpha grid cannot be empty")
    all(value -> isfinite(value) && 0 < value < 1, alpha_grid) ||
        _isoconversional_error(method, "alpha values must be finite and lie in (0, 1)")
    all(>(0), diff(alpha_grid)) ||
        _isoconversional_error(method, "alpha grid must be strictly increasing")
    if method == :advanced_vyazovkin
        minimum(alpha_grid) > configuration.advanced_vyazovkin.delta_alpha ||
            _isoconversional_error(
                method, "every alpha must exceed the configured delta_alpha"
            )
    end
    first_initial = experiments[1].initial_temperature_K
    first_final = experiments[1].final_temperature_K
    all(
        experiment ->
            experiment.initial_temperature_K == first_initial &&
            experiment.final_temperature_K == first_final,
        experiments,
    ) || _isoconversional_error(
        method, "all experiments must use the same conversion-reference temperatures"
    )
    first_composition = experiments[1].composition
    all(experiment -> experiment.composition == first_composition, experiments) ||
        _isoconversional_error(
            method, "all experiments in one analysis must have the same composition"
        )
    return nothing
end

function _linear_regression(x::Vector{Float64}, y::Vector{Float64})
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal lengths"))
    length(x) >= 3 || throw(ArgumentError("at least three observations are required"))
    centered_x = x .- mean(x)
    centered_y = y .- mean(y)
    sxx = dot(centered_x, centered_x)
    sxx > 0 || throw(ArgumentError("regression predictor has zero variance"))
    slope = dot(centered_x, centered_y) / sxx
    intercept = mean(y) - slope * mean(x)
    fitted = @. intercept + slope * x
    residuals = y .- fitted
    degrees_of_freedom = length(x) - 2
    residual_sum_squares = max(dot(residuals, residuals), 0.0)
    residual_standard_error = sqrt(residual_sum_squares / degrees_of_freedom)
    slope_standard_error = residual_standard_error / sqrt(sxx)
    total_sum_squares = dot(centered_y, centered_y)
    r_squared = if total_sum_squares <= eps(Float64)
        residual_sum_squares <= eps(Float64) ? 1.0 : NaN
    else
        1 - residual_sum_squares / total_sum_squares
    end
    return LinearRegressionDiagnostics(
        slope,
        intercept,
        slope_standard_error,
        residual_standard_error,
        r_squared,
        degrees_of_freedom,
        copy(x),
        copy(y),
        fitted,
        residuals,
    )
end

function _linear_method_variables(
    method::Symbol,
    samples::Vector{ConversionSample},
    experiments::AbstractVector{<:ProcessedExperiment},
    configuration::IsoconversionalConfig,
)
    temperatures = getproperty.(samples, :temperature_K)
    global_rates = [
        experiment.diagnostics.measured_heating_rate_K_per_min for experiment in experiments
    ]
    x = 1.0 ./ temperatures
    if method == :friedman
        rates = getproperty.(samples, :dalpha_dt_min_inv)
        all(
            rate -> isfinite(rate) && rate > configuration.friedman_minimum_rate_per_min,
            rates,
        ) || _isoconversional_error(
            method,
            "dalpha/dt must exceed $(configuration.friedman_minimum_rate_per_min) min^-1 in every experiment",
        )
        return x, log.(rates), 1.0
    elseif method == :kas
        return x, log.(global_rates ./ temperatures .^ 2), 1.0
    elseif method == :fwo
        return x, log.(global_rates), 1.052
    elseif method == :starink
        return x, log.(global_rates ./ temperatures .^ 1.92), 1.0008
    end
    return _isoconversional_error(method, "method is not a supported linear method")
end

function _linear_point(
    method::Symbol,
    alpha::Float64,
    experiments::AbstractVector{<:ProcessedExperiment},
    configuration::IsoconversionalConfig,
)
    samples = [
        interpolate_at_conversion(
            experiment, alpha; policy=configuration.interpolation_policy
        ) for experiment in experiments
    ]
    x, y, coefficient = _linear_method_variables(
        method, samples, experiments, configuration
    )
    regression = _linear_regression(x, y)
    energy = -regression.slope * GAS_CONSTANT_KJ_PER_MOL_K / coefficient
    critical = quantile(
        TDist(regression.degrees_of_freedom), 0.5 + configuration.confidence_level / 2
    )
    slope_lower = regression.slope - critical * regression.slope_standard_error
    slope_upper = regression.slope + critical * regression.slope_standard_error
    transformed = -GAS_CONSTANT_KJ_PER_MOL_K .* [slope_lower, slope_upper] ./ coefficient
    lower, upper = extrema(transformed)
    warnings = String[]
    any(sample -> sample.crossing_count > 1, samples) && push!(
        warnings,
        "one or more experiments have multiple upward conversion crossings; the earliest was used",
    )
    energy <= 0 && push!(warnings, "estimated activation energy is non-positive")
    (
        !isfinite(regression.r_squared) ||
        regression.r_squared < configuration.minimum_r_squared_warning
    ) && push!(
        warnings,
        "regression R²=$(round(regression.r_squared; digits=4)) is below the configured warning threshold $(configuration.minimum_r_squared_warning)",
    )
    status = isempty(warnings) ? :ok : :warning
    diagnostics = IsoconversionalPointDiagnostics(
        alpha,
        status,
        getproperty.(samples, :temperature_K),
        [
            experiment.diagnostics.measured_heating_rate_K_per_min for
            experiment in experiments
        ],
        getproperty.(samples, :dalpha_dt_min_inv),
        getproperty.(samples, :crossing_count),
        regression,
        nothing,
        warnings,
    )
    return energy, lower, upper, diagnostics
end

function _advanced_interval(
    experiment::ProcessedExperiment, alpha::Float64, configuration::IsoconversionalConfig
)
    advanced = configuration.advanced_vyazovkin
    alpha_nodes = collect(
        range(alpha - advanced.delta_alpha, alpha; length=advanced.integration_points)
    )
    samples = [
        interpolate_at_conversion(
            experiment, node; policy=configuration.interpolation_policy
        ) for node in alpha_nodes
    ]
    time = getproperty.(samples, :time_min)
    all(>(0), diff(time)) || _isoconversional_error(
        :advanced_vyazovkin,
        "$(experiment.source_id) does not have strictly increasing crossing times over the conversion interval ending at alpha=$alpha",
    )
    return time,
    getproperty.(samples, :temperature_K),
    maximum(getproperty.(samples, :crossing_count))
end

function _log_arrhenius_integral(
    energy_kJ_per_mol::Float64, time_min::Vector{Float64}, temperature_K::Vector{Float64}
)
    exponent = @. -energy_kJ_per_mol / (GAS_CONSTANT_KJ_PER_MOL_K * temperature_K)
    maximum_exponent = maximum(exponent)
    scaled_integral = 0.0
    for index in 2:length(time_min)
        delta_time = time_min[index] - time_min[index - 1]
        delta_time > 0 || return -Inf
        scaled_integral +=
            delta_time * (
                exp(exponent[index - 1] - maximum_exponent) +
                exp(exponent[index] - maximum_exponent)
            ) / 2
    end
    scaled_integral > 0 && isfinite(scaled_integral) || return -Inf
    return maximum_exponent + log(scaled_integral)
end

"""
    pairwise_vyazovkin_objective(log_integrals)

Evaluate `sum_i sum_(j != i) J_i/J_j` from logarithms of positive Arrhenius integrals. Equal
integrals return exactly `n*(n-1)` up to floating-point roundoff; no experiment count is
hard-coded.
"""
function pairwise_vyazovkin_objective(log_integrals::AbstractVector{<:Real})
    values = sort!(Float64.(log_integrals))
    length(values) >= 2 || return Inf
    all(isfinite, values) || return Inf
    total = 0.0
    for i in eachindex(values), j in eachindex(values)
        i == j && continue
        difference = values[i] - values[j]
        difference > log(floatmax(Float64)) && return Inf
        total += exp(difference)
    end
    return total
end

"""
    pairwise_vyazovkin_variance(log_integrals)

Evaluate the Vyazovkin–Wight pairwise ratio variance
`sum((J_i/J_j - 1)^2)/(n*(n-1))` used only for the Fisher confidence interval.
"""
function pairwise_vyazovkin_variance(log_integrals::AbstractVector{<:Real})
    values = sort!(Float64.(log_integrals))
    all(isfinite, values) || return Inf
    count = length(values) * (length(values) - 1)
    count > 0 || return Inf
    total = 0.0
    for i in eachindex(values), j in eachindex(values)
        i == j && continue
        difference = values[i] - values[j]
        difference > log(floatmax(Float64)) / 2 && return Inf
        ratio_error = exp(difference) - 1
        total += ratio_error^2
    end
    return total / count
end

function _golden_section_minimize(function_value, lower, upper, tolerance, max_iterations)
    golden_ratio = (sqrt(5.0) - 1) / 2
    left = Float64(lower)
    right = Float64(upper)
    c = right - golden_ratio * (right - left)
    d = left + golden_ratio * (right - left)
    fc = function_value(c)
    fd = function_value(d)
    evaluations = 2
    iterations = 0
    while right - left > tolerance && iterations < max_iterations
        if fc <= fd
            right = d
            d = c
            fd = fc
            c = right - golden_ratio * (right - left)
            fc = function_value(c)
        else
            left = c
            c = d
            fc = fd
            d = left + golden_ratio * (right - left)
            fd = function_value(d)
        end
        evaluations += 1
        iterations += 1
    end
    candidates = [left, c, d, right]
    candidate_values = function_value.(candidates)
    evaluations += length(candidates)
    index = argmin(candidate_values)
    return candidates[index], candidate_values[index], evaluations
end

function _scan_and_minimize(function_value, advanced::AdvancedVyazovkinConfig)
    lower, upper = advanced.energy_bounds_kJ_per_mol
    energies = collect(range(lower, upper; length=advanced.optimizer_grid_points))
    values = function_value.(energies)
    all(isfinite, values) || _isoconversional_error(
        :advanced_vyazovkin, "objective scan contains non-finite values"
    )
    best_index = argmin(values)
    bracket_lower = energies[max(1, best_index - 1)]
    bracket_upper = energies[min(length(energies), best_index + 1)]
    energy, value, evaluations = _golden_section_minimize(
        function_value,
        bracket_lower,
        bracket_upper,
        advanced.optimizer_tolerance_kJ_per_mol,
        advanced.optimizer_max_iterations,
    )
    local_minima = count(2:(length(values) - 1)) do index
        return values[index] <= values[index - 1] && values[index] <= values[index + 1]
    end
    best_index == 1 && (local_minima += 1)
    best_index == length(values) && (local_minima += 1)
    relative_range = (maximum(values) - minimum(values)) / max(abs(minimum(values)), eps())
    return (
        energy=energy,
        value=value,
        evaluations=length(values) + evaluations,
        local_minima=local_minima,
        relative_range=relative_range,
        energies=energies,
        values=values,
    )
end

function _bisect_root(function_value, lower, upper, tolerance, max_iterations)
    left = Float64(lower)
    right = Float64(upper)
    fleft = function_value(left)
    fright = function_value(right)
    fleft == 0 && return left
    fright == 0 && return right
    signbit(fleft) == signbit(fright) && return nothing
    for _ in 1:max_iterations
        midpoint = (left + right) / 2
        fmidpoint = function_value(midpoint)
        if abs(fmidpoint) <= eps(Float64) || right - left <= tolerance
            return midpoint
        end
        if signbit(fleft) == signbit(fmidpoint)
            left = midpoint
            fleft = fmidpoint
        else
            right = midpoint
        end
    end
    return (left + right) / 2
end

function _nearest_threshold_root(
    function_value,
    center::Float64,
    bound::Float64,
    direction::Symbol,
    advanced::AdvancedVyazovkinConfig,
)
    points = if direction == :lower
        collect(range(center, bound; length=advanced.optimizer_grid_points))
    else
        collect(range(center, bound; length=advanced.optimizer_grid_points))
    end
    previous_point = first(points)
    previous_value = function_value(previous_point)
    for point in Iterators.drop(points, 1)
        value = function_value(point)
        if signbit(value) != signbit(previous_value) || value == 0
            lower, upper = extrema((previous_point, point))
            return _bisect_root(
                function_value,
                lower,
                upper,
                advanced.optimizer_tolerance_kJ_per_mol,
                advanced.optimizer_max_iterations,
            ),
            false
        end
        previous_point = point
        previous_value = value
    end
    return bound, true
end

function _advanced_point(
    alpha::Float64,
    experiments::AbstractVector{<:ProcessedExperiment},
    configuration::IsoconversionalConfig,
)
    intervals = [
        _advanced_interval(experiment, alpha, configuration) for experiment in experiments
    ]
    crossing_counts = [interval[3] for interval in intervals]
    log_integrals(energy) = [
        _log_arrhenius_integral(energy, interval[1], interval[2]) for interval in intervals
    ]
    objective(energy) = pairwise_vyazovkin_objective(log_integrals(energy))
    variance(energy) = pairwise_vyazovkin_variance(log_integrals(energy))
    advanced = configuration.advanced_vyazovkin
    objective_minimum = _scan_and_minimize(objective, advanced)
    variance_minimum = _scan_and_minimize(variance, advanced)
    energy = objective_minimum.energy
    lower_bound, upper_bound = advanced.energy_bounds_kJ_per_mol
    boundary_solution =
        energy - lower_bound <= advanced.boundary_tolerance_kJ_per_mol ||
        upper_bound - energy <= advanced.boundary_tolerance_kJ_per_mol
    flat_objective = objective_minimum.relative_range <= advanced.flat_relative_tolerance

    fisher_critical = quantile(
        FDist(length(experiments) - 1, length(experiments) - 1),
        configuration.confidence_level,
    )
    fisher_variance = max(variance_minimum.value, eps(Float64))
    threshold_function(candidate_energy) =
        variance(candidate_energy) - fisher_critical * fisher_variance
    confidence_lower, lower_truncated = _nearest_threshold_root(
        threshold_function, variance_minimum.energy, lower_bound, :lower, advanced
    )
    confidence_upper, upper_truncated = _nearest_threshold_root(
        threshold_function, variance_minimum.energy, upper_bound, :upper, advanced
    )

    warnings = String[]
    any(>(1), crossing_counts) && push!(
        warnings,
        "one or more experiments have multiple upward conversion crossings; the earliest was used",
    )
    boundary_solution && push!(warnings, "objective minimum is at an energy boundary")
    flat_objective &&
        push!(warnings, "objective is flat across the configured energy range")
    objective_minimum.local_minima > 1 && push!(
        warnings,
        "objective scan contains $(objective_minimum.local_minima) local minima",
    )
    lower_truncated &&
        push!(warnings, "Fisher interval is truncated at the lower energy bound")
    upper_truncated &&
        push!(warnings, "Fisher interval is truncated at the upper energy bound")
    !(confidence_lower <= energy <= confidence_upper) && push!(
        warnings,
        "pairwise-objective estimate lies outside the variance-centered Fisher interval",
    )
    status = if boundary_solution
        :boundary_solution
    elseif flat_objective
        :flat_objective
    elseif objective_minimum.local_minima > 1
        :multimodal_objective
    elseif isempty(warnings)
        :ok
    else
        :warning
    end
    temperatures = [
        interpolate_at_conversion(
            experiment, alpha; policy=configuration.interpolation_policy
        ).temperature_K for experiment in experiments
    ]
    rates = [
        interpolate_at_conversion(
            experiment, alpha; policy=configuration.interpolation_policy
        ).dalpha_dt_min_inv for experiment in experiments
    ]
    diagnostics = AdvancedVyazovkinDiagnostics(
        objective_minimum.value,
        length(experiments) * (length(experiments) - 1),
        objective_minimum.evaluations,
        objective_minimum.local_minima,
        objective_minimum.relative_range,
        boundary_solution,
        flat_objective,
        variance_minimum.value,
        variance_minimum.energy,
        fisher_critical,
        lower_truncated,
        upper_truncated,
    )
    point = IsoconversionalPointDiagnostics(
        alpha,
        status,
        temperatures,
        [
            experiment.diagnostics.measured_heating_rate_K_per_min for
            experiment in experiments
        ],
        rates,
        crossing_counts,
        nothing,
        diagnostics,
        warnings,
    )
    return energy, confidence_lower, confidence_upper, point
end

function _analysis_fingerprint(
    method::Symbol,
    experiments::AbstractVector{<:ProcessedExperiment},
    alpha_grid::Vector{Float64},
    configuration::IsoconversionalConfig,
)
    advanced = configuration.advanced_vyazovkin
    fields = (
        method,
        alpha_grid,
        configuration.interpolation_policy,
        configuration.confidence_level,
        configuration.minimum_experiments,
        configuration.friedman_minimum_rate_per_min,
        configuration.minimum_r_squared_warning,
        advanced.delta_alpha,
        advanced.integration_points,
        advanced.energy_bounds_kJ_per_mol,
        advanced.optimizer_grid_points,
        advanced.optimizer_tolerance_kJ_per_mol,
        advanced.optimizer_max_iterations,
        advanced.boundary_tolerance_kJ_per_mol,
        advanced.flat_relative_tolerance,
        [
            (experiment.source_id, experiment.source_sha256, experiment.config_fingerprint)
            for experiment in experiments
        ],
    )
    return bytes2hex(sha256(join(repr.(fields), "|")))
end

function _failed_point(alpha::Float64, error::IsoconversionalError)
    return IsoconversionalPointDiagnostics(
        alpha,
        :invalid,
        Float64[],
        Float64[],
        Float64[],
        Int[],
        nothing,
        nothing,
        [error.message],
    )
end

"""
    analyze(method, experiments, alpha_grid; configuration)
    analyze(method, experiments, analysis_config)

Estimate conversion-dependent activation energy with `:friedman`, `:kas`, `:fwo`,
`:starink`, or `:advanced_vyazovkin`. Every experiment is interpolated at the same conversion;
unsupported points are retained as missing values with diagnostics. Linear-method intervals
are two-sided Student-t intervals for the OLS slope. Advanced-Vyazovkin intervals use the
Vyazovkin–Wight pairwise-variance F procedure and remain separate from objective quality.
"""
function analyze(
    method::Symbol,
    experiments::AbstractVector{<:ProcessedExperiment},
    alpha_grid;
    configuration::IsoconversionalConfig,
)
    grid = Float64.(collect(alpha_grid))
    _validate_analysis_inputs(method, experiments, grid, configuration)
    ordered_experiments = sort(collect(experiments); by=experiment -> experiment.source_id)
    energies = Union{Missing,Float64}[]
    lower = Union{Missing,Float64}[]
    upper = Union{Missing,Float64}[]
    point_diagnostics = IsoconversionalPointDiagnostics[]
    for alpha in grid
        try
            energy, point_lower, point_upper, diagnostics = if method == :advanced_vyazovkin
                _advanced_point(alpha, ordered_experiments, configuration)
            else
                _linear_point(method, alpha, ordered_experiments, configuration)
            end
            push!(energies, energy)
            push!(lower, point_lower)
            push!(upper, point_upper)
            push!(point_diagnostics, diagnostics)
        catch error
            if !(error isa IsoconversionalError)
                error isa ArgumentError || error isa DomainError || rethrow()
                error = IsoconversionalError(method, sprint(showerror, error))
            end
            push!(energies, missing)
            push!(lower, missing)
            push!(upper, missing)
            push!(point_diagnostics, _failed_point(alpha, error))
        end
    end
    interval_method = if method == :advanced_vyazovkin
        "Vyazovkin–Wight pairwise-variance Fisher interval"
    else
        "two-sided Student-t interval for transformed OLS slope"
    end
    return IsoconversionalResult(
        method,
        grid,
        energies,
        lower,
        upper,
        configuration.confidence_level,
        interval_method,
        [experiment.source_id for experiment in ordered_experiments],
        [experiment.config_fingerprint for experiment in ordered_experiments],
        configuration,
        point_diagnostics,
        _analysis_fingerprint(method, ordered_experiments, grid, configuration),
    )
end

function analyze(
    method::Symbol,
    experiments::AbstractVector{<:ProcessedExperiment},
    analysis_config::AnalysisConfig,
)
    return analyze(
        method,
        experiments,
        analysis_config.isoconversional.alpha_grid;
        configuration=analysis_config.isoconversional,
    )
end
