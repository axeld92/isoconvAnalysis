"""
Raised when a deconvolution request violates its data, model, or numerical contract.
"""
struct DeconvolutionError <: Exception
    experiment_id::String
    message::String
end

function Base.showerror(io::IO, error::DeconvolutionError)
    return print(io, "deconvolution failed for '$(error.experiment_id)': $(error.message)")
end

"""
Linearized uncertainty for one ordered Fraser–Suzuki peak. Tuple order is height, skew,
center, and width. Confidence limits may cross physical bounds and are diagnostics rather
than refitted profile limits.
"""
struct PeakParameterUncertainty
    standard_error::NTuple{4,Float64}
    confidence_lower::NTuple{4,Float64}
    confidence_upper::NTuple{4,Float64}
end

"""
Optimizer, reconstruction, information-criterion, area, and local-identifiability evidence
for one deconvolution fit.
"""
struct DeconvolutionDiagnostics
    status::Symbol
    converged::Bool
    optimizer_iterations::Int
    objective_evaluations::Int
    multistart_objectives::Vector{Float64}
    residual_sum_squares::Float64
    root_mean_square_error::Float64
    r_squared::Float64
    aicc::Float64
    bic::Float64
    durbin_watson::Float64
    observed_area::Float64
    reconstructed_area::Float64
    relative_area_error::Float64
    jacobian_rank::Int
    parameter_count::Int
    jacobian_condition::Float64
    maximum_absolute_parameter_correlation::Float64
    boundary_parameters::Vector{String}
    warnings::Vector{String}
end

"""
One ordered Fraser–Suzuki decomposition evaluated on the exact selected conversion interval.
Component curves and areas retain physical `K^-1` and dimensionless units respectively.
"""
struct DeconvolutionResult
    experiment_id::String
    peak_count::Int
    temperature_K::Vector{Float64}
    observed_rate_K_inv::Vector{Float64}
    fitted_rate_K_inv::Vector{Float64}
    residual_rate_K_inv::Vector{Float64}
    peaks::Vector{FraserSuzukiPeak}
    component_rate_K_inv::Vector{Vector{Float64}}
    component_areas::Vector{Float64}
    component_area_fractions::Vector{Float64}
    parameter_uncertainty::Vector{PeakParameterUncertainty}
    parameter_covariance::Matrix{Float64}
    parameter_correlation::Matrix{Float64}
    configuration::DeconvolutionConfig
    diagnostics::DeconvolutionDiagnostics
    analysis_fingerprint::String
end

"""
Complete peak-count comparison for one experiment. Every candidate result is retained,
including the raw information-criterion minimum and structural-eligibility decision.
"""
struct PeakCountComparison
    experiment_id::String
    criterion::Symbol
    criterion_minimum_peak_count::Int
    selected_peak_count::Int
    criterion_values::Vector{Float64}
    criterion_deltas::Vector{Float64}
    structurally_eligible::Vector{Bool}
    results::Vector{DeconvolutionResult}
    status::Symbol
    warnings::Vector{String}
end

"""
Diagnostics for a joint multi-rate fit with shared component skews and experiment-specific
heights, ordered centers, and widths.
"""
struct JointDeconvolutionDiagnostics
    status::Symbol
    converged::Bool
    optimizer_iterations::Int
    objective_evaluations::Int
    multistart_objectives::Vector{Float64}
    normalized_objective::Float64
    independent_normalized_objective::Float64
    relative_objective_increase::Float64
    warnings::Vector{String}
end

"""
Joint multi-rate Fraser–Suzuki result. Experiments are canonically ordered by identifier;
component skew is shared while other parameters remain experiment-specific.
"""
struct JointDeconvolutionResult
    peak_count::Int
    experiment_ids::Vector{String}
    shared_skews::Vector{Float64}
    experiment_results::Vector{DeconvolutionResult}
    configuration::DeconvolutionConfig
    diagnostics::JointDeconvolutionDiagnostics
    analysis_fingerprint::String
end

function _deconvolution_error(experiment_id, message)
    return throw(DeconvolutionError(String(experiment_id), String(message)))
end

_logistic(value) = inv(one(value) + exp(-value))

function _logit(value::Real)
    bounded = clamp(Float64(value), 1.0e-8, 1 - 1.0e-8)
    return log(bounded / (1 - bounded))
end

function _softmax_with_reference(logits)
    maximum_value = max(maximum(logits), zero(eltype(logits)))
    weights = [exp(value - maximum_value) for value in logits]
    push!(weights, exp(-maximum_value))
    return weights ./ sum(weights)
end

function _decode_peak_tuples(
    transformed,
    peak_count::Int,
    minimum_temperature::Float64,
    maximum_temperature::Float64,
    maximum_height::Float64,
    config::DeconvolutionConfig,
)
    length(transformed) == 4 * peak_count ||
        throw(DimensionMismatch("transformed peak vector has the wrong length"))
    available =
        maximum_temperature - minimum_temperature -
        (peak_count - 1) * config.minimum_center_separation_K
    available > 0 ||
        throw(ArgumentError("temperature interval is too narrow for peak count"))
    center_logits = transformed[(3 * peak_count + 1):(4 * peak_count)]
    gap_fractions = _softmax_with_reference(center_logits)
    centers = similar(center_logits, peak_count)
    cumulative_gap = zero(eltype(transformed))
    for index in 1:peak_count
        cumulative_gap += gap_fractions[index] * available
        centers[index] =
            minimum_temperature +
            cumulative_gap +
            (index - 1) * config.minimum_center_separation_K
    end
    skew_midpoint = sum(config.skew_bounds) / 2
    skew_half_range = (config.skew_bounds[2] - config.skew_bounds[1]) / 2
    width_minimum, width_maximum = config.width_bounds_K
    return [
        (
            maximum_height * _logistic(transformed[3 * index - 2]),
            skew_midpoint + skew_half_range * tanh(transformed[3 * index - 1]),
            centers[index],
            width_minimum +
            (width_maximum - width_minimum) * _logistic(transformed[3 * index]),
        ) for index in 1:peak_count
    ]
end

function _tuple_mixture(temperature, peaks)
    return [
        sum(
            _fraser_suzuki_value(value, peak[1], peak[2], peak[3], peak[4]) for
            peak in peaks
        ) for value in temperature
    ]
end

function _public_peaks(peaks)
    return [FraserSuzukiPeak(Float64.(peak)...) for peak in peaks]
end

function _encode_peaks(
    peaks::Vector{FraserSuzukiPeak},
    minimum_temperature::Float64,
    maximum_temperature::Float64,
    maximum_height::Float64,
    config::DeconvolutionConfig,
)
    peak_count = length(peaks)
    transformed = zeros(Float64, 4 * peak_count)
    skew_midpoint = sum(config.skew_bounds) / 2
    skew_half_range = (config.skew_bounds[2] - config.skew_bounds[1]) / 2
    width_minimum, width_maximum = config.width_bounds_K
    for (index, peak) in enumerate(peaks)
        transformed[3 * index - 2] = _logit(peak.height_K_inv / maximum_height)
        normalized_skew = clamp(
            (peak.skew - skew_midpoint) / skew_half_range, -1 + 1.0e-8, 1 - 1.0e-8
        )
        transformed[3 * index - 1] = atanh(normalized_skew)
        transformed[3 * index] = _logit(
            (peak.width_K - width_minimum) / (width_maximum - width_minimum)
        )
    end
    available =
        maximum_temperature - minimum_temperature -
        (peak_count - 1) * config.minimum_center_separation_K
    gaps = Float64[]
    previous_center = minimum_temperature
    for (index, peak) in enumerate(peaks)
        gap = peak.center_K - previous_center
        index > 1 && (gap -= config.minimum_center_separation_K)
        push!(gaps, max(gap, 1.0e-6 * available))
        previous_center = peak.center_K
    end
    push!(gaps, max(maximum_temperature - last(peaks).center_K, 1.0e-6 * available))
    gaps ./= sum(gaps)
    reference = gaps[end]
    transformed[(3 * peak_count + 1):(4 * peak_count)] .= log.(
        gaps[1:(end - 1)] ./ reference
    )
    return transformed
end

function _selected_curve(experiment::ProcessedExperiment, config::DeconvolutionConfig)
    lower_alpha, upper_alpha = config.fit_conversion_range
    lower_indices = _crossing_indices(experiment.analysis_alpha, lower_alpha)
    upper_indices = _crossing_indices(experiment.analysis_alpha, upper_alpha)
    isempty(lower_indices) && _deconvolution_error(
        experiment.source_id, "lower fit conversion $lower_alpha has no upward crossing"
    )
    isempty(upper_indices) && _deconvolution_error(
        experiment.source_id, "upper fit conversion $upper_alpha has no upward crossing"
    )
    lower_index = first(lower_indices)
    candidates = filter(index -> index >= lower_index, upper_indices)
    isempty(candidates) && _deconvolution_error(
        experiment.source_id,
        "upper conversion crossing precedes the selected lower crossing",
    )
    upper_index = first(candidates)
    lower_sample = interpolate_at_conversion(experiment, lower_alpha)
    upper_sample = interpolate_at_conversion(experiment, upper_alpha)
    interior = (lower_index + 1):upper_index
    temperature = vcat(
        lower_sample.temperature_K,
        experiment.temperature_K[interior],
        upper_sample.temperature_K,
    )
    rate = vcat(
        lower_sample.dalpha_dT_K_inv,
        experiment.dalpha_dT_K_inv[interior],
        upper_sample.dalpha_dT_K_inv,
    )
    all(>(0), diff(temperature)) || _deconvolution_error(
        experiment.source_id, "selected temperatures are not strictly increasing"
    )
    return temperature, rate
end

function _fit_indices(point_count::Int, maximum_fit_points::Int)
    point_count <= maximum_fit_points && return collect(1:point_count)
    return unique(round.(Int, range(1, point_count; length=maximum_fit_points)))
end

function _initial_peaks(
    temperature::Vector{Float64},
    rate::Vector{Float64},
    peak_count::Int,
    start_index::Int,
    maximum_height::Float64,
    config::DeconvolutionConfig,
)
    minimum_temperature, maximum_temperature = extrema(temperature)
    span = maximum_temperature - minimum_temperature
    centers = collect(
        range(
            minimum_temperature + span / (peak_count + 1),
            maximum_temperature - span / (peak_count + 1);
            length=peak_count,
        ),
    )
    local_maxima = [
        index for index in 2:(length(rate) - 1) if
        rate[index] >= rate[index - 1] && rate[index] >= rate[index + 1] && rate[index] > 0
    ]
    ranked = sort(local_maxima; by=index -> rate[index], rev=true)
    selected = Int[]
    for index in ranked
        all(
            abs(temperature[index] - temperature[other]) >=
            config.minimum_center_separation_K for other in selected
        ) || continue
        push!(selected, index)
        length(selected) == peak_count && break
    end
    if start_index == 1 && length(selected) == peak_count
        centers = sort(temperature[selected])
    elseif start_index > 2
        shift = 0.035 * span * sin(start_index)
        centers = clamp.(centers .+ shift, minimum_temperature + 1, maximum_temperature - 1)
    end
    base_width = clamp(
        span / (1.6 * peak_count),
        config.width_bounds_K[1] * 1.2,
        config.width_bounds_K[2] * 0.8,
    )
    peaks = FraserSuzukiPeak[]
    for (index, center) in enumerate(centers)
        nearest = argmin(abs.(temperature .- center))
        height = clamp(
            max(rate[nearest], maximum(rate) / (2 * peak_count)),
            1.0e-8,
            0.8 * maximum_height,
        )
        skew = start_index <= 2 ? 0.0 : 0.25 * sin(index + 2 * start_index)
        width = clamp(
            base_width * (1 + 0.12 * cos(index + start_index)),
            config.width_bounds_K[1] * 1.05,
            config.width_bounds_K[2] * 0.95,
        )
        push!(peaks, FraserSuzukiPeak(height, skew, center, width))
    end
    return peaks
end

function _optimize_peaks(
    temperature::Vector{Float64},
    rate::Vector{Float64},
    peak_count::Int,
    config::DeconvolutionConfig,
)
    fit_indices = _fit_indices(length(temperature), config.maximum_fit_points)
    fit_temperature = temperature[fit_indices]
    fit_rate = rate[fit_indices]
    minimum_temperature, maximum_temperature = extrema(temperature)
    rate_scale = max(maximum(abs, fit_rate), 1.0e-8)
    maximum_height = config.maximum_height_factor * rate_scale
    objective(transformed) = begin
        peaks = _decode_peak_tuples(
            transformed,
            peak_count,
            minimum_temperature,
            maximum_temperature,
            maximum_height,
            config,
        )
        residuals = _tuple_mixture(fit_temperature, peaks) .- fit_rate
        return sum(abs2, residuals) / (length(residuals) * rate_scale^2)
    end
    options = Optim.Options(;
        iterations=config.optimizer_max_iterations,
        g_tol=config.optimizer_gradient_tolerance,
        f_reltol=config.optimizer_function_tolerance,
        allow_f_increases=true,
        show_trace=false,
    )
    runs = Any[]
    objectives = Float64[]
    for start_index in 1:config.multistart_count
        initial_peaks = _initial_peaks(
            fit_temperature, fit_rate, peak_count, start_index, maximum_height, config
        )
        initial = _encode_peaks(
            initial_peaks, minimum_temperature, maximum_temperature, maximum_height, config
        )
        result = try
            Optim.optimize(
                objective,
                initial,
                Optim.LBFGS(),
                options;
                autodiff=Optim.ADTypes.AutoFiniteDiff(),
            )
        catch
            nothing
        end
        isnothing(result) && continue
        value = Optim.minimum(result)
        isfinite(value) || continue
        push!(runs, result)
        push!(objectives, value)
    end
    isempty(runs) && _deconvolution_error("optimizer", "every multistart run failed")
    best_index = argmin(objectives)
    best = runs[best_index]
    tuples = _decode_peak_tuples(
        Optim.minimizer(best),
        peak_count,
        minimum_temperature,
        maximum_temperature,
        maximum_height,
        config,
    )
    return (
        peaks=_public_peaks(tuples),
        result=best,
        objectives=objectives,
        maximum_height=maximum_height,
    )
end

function _physical_parameter_vector(peaks::Vector{FraserSuzukiPeak})
    parameters = Float64[]
    for peak in peaks
        append!(parameters, (peak.height_K_inv, peak.skew, peak.center_K, peak.width_K))
    end
    return parameters
end

function _peaks_from_parameter_vector(parameters::Vector{Float64})
    length(parameters) % 4 == 0 || throw(DimensionMismatch("invalid peak parameter vector"))
    return [
        FraserSuzukiPeak(parameters[index:(index + 3)]...) for
        index in 1:4:length(parameters)
    ]
end

function _model_jacobian(temperature, peaks)
    parameters = _physical_parameter_vector(peaks)
    jacobian = zeros(Float64, length(temperature), length(parameters))
    for index in eachindex(parameters)
        step =
            cbrt(eps(Float64)) * max(abs(parameters[index]), index % 4 == 1 ? 1.0e-3 : 1.0)
        upper = copy(parameters)
        lower = copy(parameters)
        upper[index] += step
        lower[index] -= step
        if index % 4 == 1 && lower[index] < 0
            lower[index] = parameters[index]
            denominator = step
        elseif index % 4 == 0 && lower[index] <= 0
            lower[index] = parameters[index]
            denominator = step
        else
            denominator = 2 * step
        end
        upper_curve = fraser_suzuki_mixture(
            temperature, _peaks_from_parameter_vector(upper)
        )
        lower_curve = fraser_suzuki_mixture(
            temperature, _peaks_from_parameter_vector(lower)
        )
        jacobian[:, index] .= (upper_curve .- lower_curve) ./ denominator
    end
    return jacobian
end

function _parameter_uncertainty(temperature, residuals, peaks, confidence_level)
    parameters = _physical_parameter_vector(peaks)
    jacobian = _model_jacobian(temperature, peaks)
    degrees_of_freedom = length(temperature) - length(parameters)
    variance = sum(abs2, residuals) / degrees_of_freedom
    information = transpose(jacobian) * jacobian
    covariance = variance .* pinv(information; rtol=sqrt(eps(Float64)))
    standard_errors = sqrt.(max.(diag(covariance), 0.0))
    critical = quantile(TDist(degrees_of_freedom), 0.5 + confidence_level / 2)
    lower = parameters .- critical .* standard_errors
    upper = parameters .+ critical .* standard_errors
    uncertainties = [
        PeakParameterUncertainty(
            Tuple(standard_errors[index:(index + 3)]),
            Tuple(lower[index:(index + 3)]),
            Tuple(upper[index:(index + 3)]),
        ) for index in 1:4:length(parameters)
    ]
    scale = sqrt.(max.(diag(covariance), 0.0))
    correlation = zeros(Float64, size(covariance))
    for row in axes(covariance, 1), column in axes(covariance, 2)
        denominator = scale[row] * scale[column]
        correlation[row, column] =
            denominator > 0 ? covariance[row, column] / denominator : 0.0
    end
    for index in axes(correlation, 1)
        correlation[index, index] = 1.0
    end
    jacobian_rank = rank(jacobian; rtol=sqrt(eps(Float64)))
    column_norms = [norm(view(jacobian, :, index)) for index in axes(jacobian, 2)]
    normalized_jacobian = copy(jacobian)
    for index in axes(normalized_jacobian, 2)
        column_norms[index] > 0 && (normalized_jacobian[:, index] ./= column_norms[index])
    end
    singular_values = svdvals(normalized_jacobian)
    jacobian_condition = if isempty(singular_values) || last(singular_values) <= 0
        Inf
    else
        first(singular_values) / last(singular_values)
    end
    off_diagonal = [
        abs(correlation[row, column]) for
        row in axes(correlation, 1), column in axes(correlation, 2) if row != column
    ]
    maximum_correlation = isempty(off_diagonal) ? 0.0 : maximum(off_diagonal)
    return (
        uncertainty=uncertainties,
        covariance=covariance,
        correlation=correlation,
        rank=jacobian_rank,
        condition=jacobian_condition,
        maximum_correlation=maximum_correlation,
    )
end

function _boundary_parameters(
    peaks, temperature, maximum_height, config::DeconvolutionConfig
)
    names = String[]
    fraction = config.boundary_fraction_tolerance
    temperature_minimum, temperature_maximum = extrema(temperature)
    for (index, peak) in enumerate(peaks)
        peak.height_K_inv <= fraction * maximum_height &&
            push!(names, "peak_$index.height_lower")
        peak.height_K_inv >= (1 - fraction) * maximum_height &&
            push!(names, "peak_$index.height_upper")
        skew_fraction =
            (peak.skew - config.skew_bounds[1]) / diff(collect(config.skew_bounds))[1]
        skew_fraction <= fraction && push!(names, "peak_$index.skew_lower")
        skew_fraction >= 1 - fraction && push!(names, "peak_$index.skew_upper")
        width_fraction =
            (peak.width_K - config.width_bounds_K[1]) /
            diff(collect(config.width_bounds_K))[1]
        width_fraction <= fraction && push!(names, "peak_$index.width_lower")
        width_fraction >= 1 - fraction && push!(names, "peak_$index.width_upper")
        if index == 1
            peak.center_K - temperature_minimum <=
            fraction * (temperature_maximum - temperature_minimum) &&
                push!(names, "peak_$index.center_lower")
        else
            separation = peak.center_K - peaks[index - 1].center_K
            separation - config.minimum_center_separation_K <=
            fraction * (temperature_maximum - temperature_minimum) &&
                push!(names, "peak_$index.center_separation")
        end
    end
    return names
end

function _deconvolution_fingerprint(experiment, peak_count, config)
    fields = (
        experiment.source_id,
        experiment.source_sha256,
        experiment.config_fingerprint,
        peak_count,
        config.peak_counts,
        config.selection_criterion,
        config.fit_conversion_range,
        config.maximum_fit_points,
        config.minimum_center_separation_K,
        config.skew_bounds,
        config.width_bounds_K,
        config.maximum_height_factor,
        config.minimum_component_area_fraction,
        config.multistart_count,
        config.optimizer_max_iterations,
        config.optimizer_gradient_tolerance,
        config.optimizer_function_tolerance,
        config.identifiability_condition_warning,
        config.boundary_fraction_tolerance,
        config.confidence_level,
    )
    return bytes2hex(sha256(join(repr.(fields), "|")))
end

function _build_deconvolution_result(
    experiment::ProcessedExperiment,
    peak_count::Int,
    temperature,
    observed,
    optimized,
    config::DeconvolutionConfig,
    ;
    fingerprint::Union{Nothing,String}=nothing,
)
    peaks = optimized.peaks
    components = [fraser_suzuki(temperature, peak) for peak in peaks]
    fitted = reduce(+, components)
    residuals = observed .- fitted
    residual_sum_squares = sum(abs2, residuals)
    root_mean_square_error = sqrt(residual_sum_squares / length(residuals))
    centered = observed .- mean(observed)
    total_sum_squares = sum(abs2, centered)
    r_squared = total_sum_squares > 0 ? 1 - residual_sum_squares / total_sum_squares : NaN
    parameter_count = 4 * peak_count
    point_count = length(observed)
    criterion_indices = _fit_indices(point_count, config.maximum_fit_points)
    criterion_point_count = length(criterion_indices)
    criterion_rss = max(sum(abs2, residuals[criterion_indices]), eps(Float64))
    aic =
        criterion_point_count * log(criterion_rss / criterion_point_count) +
        2 * parameter_count
    aicc = if criterion_point_count > parameter_count + 1
        aic +
        2 * parameter_count * (parameter_count + 1) /
        (criterion_point_count - parameter_count - 1)
    else
        Inf
    end
    bic =
        criterion_point_count * log(criterion_rss / criterion_point_count) +
        parameter_count * log(criterion_point_count)
    difference_residual = diff(residuals)
    durbin_watson = if residual_sum_squares > 0
        sum(abs2, difference_residual) / residual_sum_squares
    else
        2.0
    end
    component_areas = [
        last(cumulative_trapezoid(temperature, curve)) for curve in components
    ]
    reconstructed_area = sum(component_areas)
    observed_area = last(cumulative_trapezoid(temperature, observed))
    area_fractions = if reconstructed_area > 0
        component_areas ./ reconstructed_area
    else
        fill(NaN, peak_count)
    end
    relative_area_error =
        abs(reconstructed_area - observed_area) / max(abs(observed_area), eps(Float64))
    uncertainty = _parameter_uncertainty(
        temperature, residuals, peaks, config.confidence_level
    )
    boundary_parameters = _boundary_parameters(
        peaks, temperature, optimized.maximum_height, config
    )
    warnings = String[]
    converged = Optim.converged(optimized.result)
    converged || push!(warnings, "best optimizer run did not satisfy convergence criteria")
    any(<(0), observed) && push!(
        warnings,
        "observed derivative contains negative values; the nonnegative peak model does not clip them",
    )
    !isempty(boundary_parameters) &&
        push!(warnings, "one or more peak parameters lie near configured bounds")
    uncertainty.rank < parameter_count &&
        push!(warnings, "model Jacobian is rank deficient")
    uncertainty.condition > config.identifiability_condition_warning && push!(
        warnings,
        "Jacobian condition $(uncertainty.condition) exceeds the identifiability warning threshold",
    )
    uncertainty.maximum_correlation > 0.98 &&
        push!(warnings, "at least one linearized parameter correlation exceeds 0.98")
    relative_area_error > 0.02 && push!(
        warnings,
        "reconstructed derivative area differs from the observed area by more than 2%",
    )
    status = if !converged
        :not_converged
    elseif uncertainty.rank < parameter_count ||
        uncertainty.condition > config.identifiability_condition_warning
        :nonidentifiable
    elseif isempty(warnings)
        :ok
    else
        :warning
    end
    diagnostics = DeconvolutionDiagnostics(
        status,
        converged,
        Optim.iterations(optimized.result),
        Optim.f_calls(optimized.result),
        sort(optimized.objectives),
        residual_sum_squares,
        root_mean_square_error,
        r_squared,
        aicc,
        bic,
        durbin_watson,
        observed_area,
        reconstructed_area,
        relative_area_error,
        uncertainty.rank,
        parameter_count,
        uncertainty.condition,
        uncertainty.maximum_correlation,
        boundary_parameters,
        warnings,
    )
    return DeconvolutionResult(
        experiment.source_id,
        peak_count,
        temperature,
        observed,
        fitted,
        residuals,
        peaks,
        components,
        component_areas,
        area_fractions,
        uncertainty.uncertainty,
        uncertainty.covariance,
        uncertainty.correlation,
        config,
        diagnostics,
        if isnothing(fingerprint)
            _deconvolution_fingerprint(experiment, peak_count, config)
        else
            fingerprint
        end,
    )
end

"""
    fit_deconvolution(experiment, peak_count; configuration)
    fit_deconvolution(experiment, peak_count, analysis_config)

Fit an ordered, nonnegative Fraser–Suzuki mixture to `dalpha/dT` over the configured exact
conversion interval. Bounds are enforced by a smooth transformation and deterministic
multistart L-BFGS optimization. The returned result always preserves physical units and
contains reconstruction, information-criterion, area, uncertainty, boundary, and
identifiability diagnostics.
"""
function fit_deconvolution(
    experiment::ProcessedExperiment, peak_count::Integer; configuration::DeconvolutionConfig
)
    count = Int(peak_count)
    count in configuration.peak_counts || _deconvolution_error(
        experiment.source_id,
        "peak count $count is not configured; allowed counts are $(configuration.peak_counts)",
    )
    temperature, observed = _selected_curve(experiment, configuration)
    length(observed) > 4 * count + 2 || _deconvolution_error(
        experiment.source_id, "too few selected points for $count peaks"
    )
    span = last(temperature) - first(temperature)
    span > (count - 1) * configuration.minimum_center_separation_K ||
        _deconvolution_error(experiment.source_id, "temperature interval is too narrow")
    optimized = _optimize_peaks(temperature, observed, count, configuration)
    return _build_deconvolution_result(
        experiment, count, temperature, observed, optimized, configuration
    )
end

function fit_deconvolution(
    experiment::ProcessedExperiment, peak_count::Integer, config::AnalysisConfig
)
    return fit_deconvolution(experiment, peak_count; configuration=config.deconvolution)
end

"""
    compare_peak_counts(experiment; configuration)

Fit every configured peak count and calculate BIC or AICc. Selection excludes candidates
with nonconvergence, local non-identifiability, active physical bounds, or a component below
the configured minimum area fraction. The raw criterion minimum and every candidate fit are
retained for review.
"""
function compare_peak_counts(
    experiment::ProcessedExperiment; configuration::DeconvolutionConfig
)
    results = [
        fit_deconvolution(experiment, count; configuration=configuration) for
        count in configuration.peak_counts
    ]
    values = if configuration.selection_criterion == :bic
        [result.diagnostics.bic for result in results]
    else
        [result.diagnostics.aicc for result in results]
    end
    criterion_index = argmin(values)
    deltas = values .- values[criterion_index]
    eligible = [
        result.diagnostics.converged &&
            result.diagnostics.jacobian_rank == result.diagnostics.parameter_count &&
            result.diagnostics.jacobian_condition <=
            configuration.identifiability_condition_warning &&
            isempty(result.diagnostics.boundary_parameters) &&
            minimum(result.component_area_fractions) >=
            configuration.minimum_component_area_fraction for result in results
    ]
    eligible_indices = findall(eligible)
    ordered_deltas = sort(deltas)
    warnings = String[]
    selected_index = if isempty(eligible_indices)
        push!(warnings, "no peak-count candidate passed structural eligibility checks")
        criterion_index
    else
        eligible_indices[argmin(values[eligible_indices])]
    end
    status = if isempty(eligible_indices)
        :warning
    elseif selected_index != criterion_index
        push!(warnings, "raw criterion minimum was rejected by structural eligibility checks")
        :constraint_filtered
    elseif length(ordered_deltas) > 1 && ordered_deltas[2] < 2
        push!(
            warnings,
            "best and second-best peak counts differ by less than two criterion units",
        )
        :ambiguous
    else
        :selected
    end
    results[selected_index].diagnostics.status in (:not_converged, :nonidentifiable) &&
        begin
            push!(
                warnings,
                "selected peak-count fit has status $(results[selected_index].diagnostics.status)",
            )
            status = :warning
        end
    return PeakCountComparison(
        experiment.source_id,
        configuration.selection_criterion,
        results[criterion_index].peak_count,
        results[selected_index].peak_count,
        values,
        deltas,
        eligible,
        results,
        status,
        warnings,
    )
end

function compare_peak_counts(experiment::ProcessedExperiment, config::AnalysisConfig)
    return compare_peak_counts(experiment; configuration=config.deconvolution)
end

function _skew_to_transformed(skew::Float64, config::DeconvolutionConfig)
    midpoint = sum(config.skew_bounds) / 2
    half_range = (config.skew_bounds[2] - config.skew_bounds[1]) / 2
    normalized = clamp((skew - midpoint) / half_range, -1 + 1.0e-8, 1 - 1.0e-8)
    return atanh(normalized)
end

function _decode_joint_peak_tuples(
    transformed,
    shared_transformed,
    peak_count,
    minimum_temperature,
    maximum_temperature,
    maximum_height,
    config,
)
    length(transformed) == 3 * peak_count ||
        throw(DimensionMismatch("joint experiment block has the wrong length"))
    independent_transformed = similar(transformed, 4 * peak_count)
    for peak_index in 1:peak_count
        independent_transformed[3 * peak_index - 2] = transformed[2 * peak_index - 1]
        independent_transformed[3 * peak_index - 1] = shared_transformed[peak_index]
        independent_transformed[3 * peak_index] = transformed[2 * peak_index]
    end
    independent_transformed[(3 * peak_count + 1):(4 * peak_count)] .= transformed[(2 * peak_count + 1):(3 * peak_count)]
    return _decode_peak_tuples(
        independent_transformed,
        peak_count,
        minimum_temperature,
        maximum_temperature,
        maximum_height,
        config,
    )
end

function _joint_initial_vector(independent_results, maximum_heights, config)
    peak_count = first(independent_results).peak_count
    experiment_count = length(independent_results)
    shared_skews = [
        median(result.peaks[index].skew for result in independent_results) for
        index in 1:peak_count
    ]
    transformed = zeros(Float64, peak_count + 3 * peak_count * experiment_count)
    transformed[1:peak_count] .= _skew_to_transformed.(shared_skews, Ref(config))
    for experiment_index in 1:experiment_count
        result = independent_results[experiment_index]
        encoded = _encode_peaks(
            result.peaks,
            first(result.temperature_K),
            last(result.temperature_K),
            maximum_heights[experiment_index],
            config,
        )
        block_start = peak_count + (experiment_index - 1) * 3 * peak_count
        for peak_index in 1:peak_count
            transformed[block_start + 2 * peak_index - 1] = encoded[3 * peak_index - 2]
            transformed[block_start + 2 * peak_index] = encoded[3 * peak_index]
        end
        transformed[(block_start + 2 * peak_count + 1):(block_start + 3 * peak_count)] .= encoded[(3 * peak_count + 1):(4 * peak_count)]
    end
    return transformed
end

function _joint_fingerprint(experiments, peak_count, config)
    fields = (
        peak_count,
        [
            (experiment.source_id, experiment.source_sha256, experiment.config_fingerprint)
            for experiment in experiments
        ],
        config.peak_counts,
        config.selection_criterion,
        config.fit_conversion_range,
        config.maximum_fit_points,
        config.minimum_center_separation_K,
        config.skew_bounds,
        config.width_bounds_K,
        config.maximum_height_factor,
        config.minimum_component_area_fraction,
        config.multistart_count,
        config.joint_multistart_count,
        config.optimizer_max_iterations,
        config.optimizer_gradient_tolerance,
        config.optimizer_function_tolerance,
        config.identifiability_condition_warning,
        config.boundary_fraction_tolerance,
        config.confidence_level,
    )
    return bytes2hex(sha256(join(repr.(fields), "|")))
end

"""
    fit_joint_deconvolution(experiments, peak_count; configuration)

Jointly refine a common ordered component decomposition across heating programs. Skew is
shared by component; height, center, and width remain experiment-specific. Independent fits
provide deterministic initial values and a nested objective benchmark. The returned
experiment order and fingerprint are invariant to caller order.
"""
function fit_joint_deconvolution(
    experiments::AbstractVector{<:ProcessedExperiment},
    peak_count::Integer;
    configuration::DeconvolutionConfig,
    independent_results::Union{Nothing,AbstractVector{<:DeconvolutionResult}}=nothing,
)
    count = Int(peak_count)
    count in configuration.peak_counts || _deconvolution_error(
        "joint",
        "peak count $count is not configured; allowed counts are $(configuration.peak_counts)",
    )
    length(experiments) >= 2 ||
        _deconvolution_error("joint", "at least two experiments are required")
    length(unique(experiment.source_id for experiment in experiments)) ==
    length(experiments) ||
        _deconvolution_error("joint", "experiment identifiers must be unique")
    first_composition = first(experiments).composition
    all(experiment -> experiment.composition == first_composition, experiments) ||
        _deconvolution_error("joint", "all experiments must have the same composition")
    ordered = sort(collect(experiments); by=experiment -> experiment.source_id)
    curves = [_selected_curve(experiment, configuration) for experiment in ordered]
    independent = if isnothing(independent_results)
        [
            fit_deconvolution(experiment, count; configuration=configuration) for
            experiment in ordered
        ]
    else
        supplied = sort(collect(independent_results); by=result -> result.experiment_id)
        [result.experiment_id for result in supplied] == [experiment.source_id for experiment in ordered] ||
            _deconvolution_error(
                "joint", "supplied independent results do not match the experiments"
            )
        all(result -> result.peak_count == count, supplied) || _deconvolution_error(
            "joint", "supplied independent results use a different peak count"
        )
        all(zip(supplied, ordered)) do (result, experiment)
            return result.analysis_fingerprint ==
                   _deconvolution_fingerprint(experiment, count, configuration)
        end || _deconvolution_error(
            "joint", "supplied independent results use a different configuration"
        )
        supplied
    end
    maximum_heights = [
        configuration.maximum_height_factor * max(maximum(abs, curve[2]), 1.0e-8) for
        curve in curves
    ]
    scales = [max(maximum(abs, curve[2]), 1.0e-8) for curve in curves]
    fit_indices = [
        _fit_indices(length(curve[1]), configuration.maximum_fit_points) for curve in curves
    ]
    initial = _joint_initial_vector(independent, maximum_heights, configuration)
    experiment_count = length(ordered)
    objective(transformed) = begin
        shared = view(transformed, 1:count)
        total = zero(eltype(transformed))
        for experiment_index in 1:experiment_count
            block_start = count + (experiment_index - 1) * 3 * count
            block = view(transformed, (block_start + 1):(block_start + 3 * count))
            temperature, observed = curves[experiment_index]
            indices = fit_indices[experiment_index]
            tuples = _decode_joint_peak_tuples(
                block,
                shared,
                count,
                first(temperature),
                last(temperature),
                maximum_heights[experiment_index],
                configuration,
            )
            residuals =
                _tuple_mixture(temperature[indices], tuples) .- observed[indices]
            total +=
                sum(abs2, residuals) / (length(indices) * scales[experiment_index]^2)
        end
        return total / experiment_count
    end
    options = Optim.Options(;
        iterations=configuration.optimizer_max_iterations,
        g_tol=configuration.optimizer_gradient_tolerance,
        f_reltol=configuration.optimizer_function_tolerance,
        allow_f_increases=true,
        show_trace=false,
    )
    runs = Any[]
    objectives = Float64[]
    for start_index in 1:configuration.joint_multistart_count
        start = copy(initial)
        if start_index > 1
            start[1:count] .+= [0.15 * sin(start_index + index) for index in 1:count]
            for index in (count + 2 * count + 1):(3 * count):length(start)
                stop = min(index + count - 1, length(start))
                start[index:stop] .+= 0.1 * cos(start_index)
            end
        end
        result = try
            Optim.optimize(
                objective,
                start,
                Optim.LBFGS(),
                options;
                autodiff=Optim.ADTypes.AutoFiniteDiff(),
            )
        catch
            nothing
        end
        isnothing(result) && continue
        value = Optim.minimum(result)
        isfinite(value) || continue
        push!(runs, result)
        push!(objectives, value)
    end
    isempty(runs) && _deconvolution_error("joint", "every joint multistart run failed")
    best = runs[argmin(objectives)]
    best_transformed = Optim.minimizer(best)
    shared_transformed = view(best_transformed, 1:count)
    skew_midpoint = sum(configuration.skew_bounds) / 2
    skew_half_range = diff(collect(configuration.skew_bounds))[1] / 2
    shared_skews = skew_midpoint .+ skew_half_range .* tanh.(shared_transformed)
    fingerprint = _joint_fingerprint(ordered, count, configuration)
    results = DeconvolutionResult[]
    for experiment_index in 1:experiment_count
        block_start = count + (experiment_index - 1) * 3 * count
        block = view(best_transformed, (block_start + 1):(block_start + 3 * count))
        temperature, observed = curves[experiment_index]
        tuples = _decode_joint_peak_tuples(
            block,
            shared_transformed,
            count,
            first(temperature),
            last(temperature),
            maximum_heights[experiment_index],
            configuration,
        )
        optimized = (
            peaks=_public_peaks(tuples),
            result=best,
            objectives=objectives,
            maximum_height=maximum_heights[experiment_index],
        )
        push!(
            results,
            _build_deconvolution_result(
                ordered[experiment_index],
                count,
                temperature,
                observed,
                optimized,
                configuration;
                fingerprint=bytes2hex(
                    sha256("$fingerprint|$(ordered[experiment_index].source_id)")
                ),
            ),
        )
    end
    independent_objective = mean(
        sum(abs2, independent[index].residual_rate_K_inv[fit_indices[index]]) /
        (length(fit_indices[index]) * scales[index]^2) for index in eachindex(independent)
    )
    normalized_objective = Optim.minimum(best)
    relative_increase =
        (normalized_objective - independent_objective) /
        max(independent_objective, eps(Float64))
    warnings = String[]
    converged = Optim.converged(best)
    converged || push!(warnings, "joint optimizer did not satisfy convergence criteria")
    relative_increase > 0.05 && push!(
        warnings, "shared-skew joint objective exceeds independent fits by more than 5%"
    )
    any(result -> result.diagnostics.status == :nonidentifiable, results) && push!(
        warnings, "one or more conditional experiment fits are locally non-identifiable"
    )
    status = if !converged
        :not_converged
    elseif isempty(warnings)
        :ok
    else
        :warning
    end
    diagnostics = JointDeconvolutionDiagnostics(
        status,
        converged,
        Optim.iterations(best),
        Optim.f_calls(best),
        sort(objectives),
        normalized_objective,
        independent_objective,
        relative_increase,
        warnings,
    )
    return JointDeconvolutionResult(
        count,
        [experiment.source_id for experiment in ordered],
        shared_skews,
        results,
        configuration,
        diagnostics,
        fingerprint,
    )
end

function fit_joint_deconvolution(
    experiments::AbstractVector{<:ProcessedExperiment},
    peak_count::Integer,
    config::AnalysisConfig,
)
    return fit_joint_deconvolution(
        experiments, peak_count; configuration=config.deconvolution
    )
end
