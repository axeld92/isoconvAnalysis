"""
Raised when a raw experiment cannot be transformed under an explicit preprocessing
configuration.
"""
struct PreprocessingError <: Exception
    experiment_id::String
    message::String
end

function Base.showerror(io::IO, error::PreprocessingError)
    return print(io, "preprocessing failed for '$(error.experiment_id)': $(error.message)")
end

"""
Auditable counts, reference values, numerical checks, and warnings for one preprocessing
result. Counts always refer to rows in the raw experiment unless their field name states
otherwise.
"""
struct PreprocessingDiagnostics
    input_row_count::Int
    finite_row_count::Int
    segment_row_count::Int
    retained_row_count::Int
    dropped_invalid_row_count::Int
    dropped_segment_row_count::Int
    dropped_interval_row_count::Int
    selected_segment::Union{Nothing,Int}
    initial_reference_point_count::Int
    final_reference_point_count::Int
    initial_mass_percent::Float64
    final_mass_percent::Float64
    measured_heating_rate_K_per_min::Float64
    alpha_minimum::Float64
    alpha_maximum::Float64
    alpha_outside_unit_interval_count::Int
    conversion_reversal_count::Int
    negative_derivative_count::Int
    reconstruction_rmse::Float64
    reconstruction_max_abs_error::Float64
    warnings::Vector{String}
end

"""
An immutable preprocessing record whose arrays are newly allocated and do not alias the raw
`Experiment`. `alpha` is the unsmoothed conversion; `analysis_alpha` is the configured curve
used for differentiation. Temperature is absolute, derivative units are encoded in field
names, and `source_row_indices` map every retained value back to its acquisition row.
"""
struct ProcessedExperiment
    source_id::String
    source_file::String
    source_variable::String
    source_sha256::String
    composition::Composition
    source_row_indices::Vector{Int}
    temperature_K::Vector{Float64}
    raw_time_min::Vector{Float64}
    time_min::Vector{Float64}
    mass_percent::Vector{Float64}
    mass_fraction::Vector{Float64}
    alpha::Vector{Float64}
    analysis_alpha::Vector{Float64}
    dalpha_dT_K_inv::Vector{Float64}
    dalpha_dt_min_inv::Vector{Float64}
    heating_rate_K_per_min::Vector{Float64}
    reconstructed_alpha::Vector{Float64}
    reconstruction_residual::Vector{Float64}
    initial_temperature_K::Float64
    final_temperature_K::Float64
    preprocessing::PreprocessingConfig
    config_fingerprint::String
    diagnostics::PreprocessingDiagnostics
end

function _preprocessing_error(experiment::Experiment, message::AbstractString)
    return throw(PreprocessingError(experiment.id, String(message)))
end

"""
    finite_difference_derivative(x, y)

Differentiate samples on a strictly increasing, possibly irregular grid. A local quadratic
through three adjacent points is used at both boundaries and in the interior, so quadratic
polynomials are reproduced to floating-point accuracy.
"""
function finite_difference_derivative(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal lengths"))
    length(x) >= 3 || throw(ArgumentError("at least three points are required"))
    x_values = Float64.(x)
    y_values = Float64.(y)
    all(isfinite, x_values) && all(isfinite, y_values) ||
        throw(ArgumentError("x and y must contain only finite values"))
    all(>(0), diff(x_values)) || throw(ArgumentError("x must be strictly increasing"))

    derivative = similar(y_values)
    last_index = length(x_values)
    for index in eachindex(x_values)
        indices = if index == 1
            1:3
        elseif index == last_index
            (last_index - 2):last_index
        else
            (index - 1):(index + 1)
        end
        centered = x_values[indices] .- x_values[index]
        design = hcat(ones(3), centered, centered .^ 2)
        derivative[index] = (design \ y_values[indices])[2]
    end
    return derivative
end

function _local_indices(
    sorted_coordinate::Vector{Float64},
    order::Vector{Int},
    coordinate::Vector{Float64},
    index::Int,
    half_window::Float64,
    minimum_points::Int,
)
    center = coordinate[index]
    lower = searchsortedfirst(sorted_coordinate, center - half_window)
    upper = searchsortedlast(sorted_coordinate, center + half_window)
    if upper >= lower && upper - lower + 1 >= minimum_points
        return order[lower:upper]
    end
    count = min(minimum_points, length(coordinate))
    nearest = partialsortperm(abs.(coordinate .- center), 1:count)
    return sort(nearest)
end

"""
    local_polynomial_estimate(x, y, neighborhood; half_window, degree=3)

Return `(fitted, derivative)` from weighted local-polynomial regressions of `y` against
strictly increasing `x`. Neighborhoods are selected in the physical coordinate supplied by
`neighborhood` (temperature in the preprocessing pipeline), and `half_window` has that
coordinate's unit. The derivative is with respect to `x`.
"""
function local_polynomial_estimate(
    x::AbstractVector,
    y::AbstractVector,
    neighborhood::AbstractVector;
    half_window::Real,
    degree::Integer=3,
)
    length(x) == length(y) == length(neighborhood) ||
        throw(DimensionMismatch("x, y, and neighborhood must have equal lengths"))
    1 <= degree <= 5 || throw(ArgumentError("degree must be between 1 and 5"))
    half_window > 0 || throw(ArgumentError("half_window must be positive"))
    minimum_points = Int(degree) + 2
    length(x) >= minimum_points ||
        throw(ArgumentError("at least $minimum_points points are required"))

    x_values = Float64.(x)
    y_values = Float64.(y)
    coordinate = Float64.(neighborhood)
    all(isfinite, x_values) && all(isfinite, y_values) && all(isfinite, coordinate) ||
        throw(ArgumentError("local-polynomial inputs must be finite"))
    all(>(0), diff(x_values)) || throw(ArgumentError("x must be strictly increasing"))

    window = Float64(half_window)
    order = sortperm(coordinate)
    sorted_coordinate = coordinate[order]
    fitted = similar(y_values)
    derivative = similar(y_values)

    for index in eachindex(x_values)
        indices = _local_indices(
            sorted_coordinate, order, coordinate, index, window, minimum_points
        )
        centered_x = x_values[indices] .- x_values[index]
        scale_x = maximum(abs, centered_x)
        scale_x > 0 || throw(ArgumentError("local window contains no distinct x values"))
        scaled_x = centered_x ./ scale_x

        raw_distance = abs.(coordinate[indices] .- coordinate[index])
        effective_window = max(window, maximum(raw_distance))
        weights = @. exp(-0.5 * (raw_distance / effective_window)^2)
        design = reduce(hcat, (scaled_x .^ power for power in 0:Int(degree)))
        square_root_weight = sqrt.(weights)
        weighted_design = design .* square_root_weight
        coefficients = weighted_design \ (y_values[indices] .* square_root_weight)
        fitted[index] = coefficients[1]
        derivative[index] = coefficients[2] / scale_x
    end
    return fitted, derivative
end

function _robust_local_linear_reference(
    temperature_K::Vector{Float64},
    mass_percent::Vector{Float64},
    reference_K::Float64,
    half_window_K::Float64,
)
    indices = findall(value -> abs(value - reference_K) <= half_window_K, temperature_K)
    length(indices) >= 3 || throw(
        ArgumentError(
            "reference at $(reference_K - 273.15) °C has fewer than three points within ±$half_window_K K",
        ),
    )
    centered = temperature_K[indices] .- reference_K
    design = hcat(ones(length(indices)), centered)
    response = mass_percent[indices]
    coefficients = design \ response

    for _ in 1:12
        residual = response - design * coefficients
        residual_center = median(residual)
        scale = 1.4826 * median(abs.(residual .- residual_center))
        scale <= eps(Float64) && break
        standardized = abs.(residual .- residual_center) ./ (1.345 * scale)
        weights = @. ifelse(standardized <= 1, 1.0, 1 / standardized)
        square_root_weight = sqrt.(weights)
        updated = (design .* square_root_weight) \ (response .* square_root_weight)
        norm(updated - coefficients) <= 1.0e-12 * max(norm(coefficients), 1.0) &&
            return updated[1], length(indices)
        coefficients = updated
    end
    return coefficients[1], length(indices)
end

function _linear_interpolation_reference(
    temperature_K::Vector{Float64}, mass_percent::Vector{Float64}, reference_K::Float64
)
    minimum_temperature, maximum_temperature = extrema(temperature_K)
    minimum_temperature <= reference_K <= maximum_temperature ||
        throw(ArgumentError("reference temperature lies outside the selected ramp segment"))

    exact = findall(
        value -> isapprox(value, reference_K; atol=1.0e-10, rtol=0), temperature_K
    )
    !isempty(exact) && return median(mass_percent[exact]), length(exact)

    below_temperature = maximum(filter(<(reference_K), temperature_K))
    above_temperature = minimum(filter(>(reference_K), temperature_K))
    below_mass = median(mass_percent[temperature_K .== below_temperature])
    above_mass = median(mass_percent[temperature_K .== above_temperature])
    fraction = (reference_K - below_temperature) / (above_temperature - below_temperature)
    return below_mass + fraction * (above_mass - below_mass), 2
end

function _reference_mass(
    experiment::Experiment,
    temperature_K::Vector{Float64},
    mass_percent::Vector{Float64},
    reference_K::Float64,
    config::PreprocessingConfig,
)
    try
        if config.reference_mass_method == :local_robust_linear
            return _robust_local_linear_reference(
                temperature_K, mass_percent, reference_K, config.reference_half_window_K
            )
        elseif config.reference_mass_method == :linear_interpolation
            return _linear_interpolation_reference(temperature_K, mass_percent, reference_K)
        end
    catch error
        error isa ArgumentError || rethrow()
        _preprocessing_error(experiment, sprint(showerror, error))
    end
    return _preprocessing_error(
        experiment, "unsupported reference-mass method $(config.reference_mass_method)"
    )
end

function _select_segment_rows(
    experiment::Experiment,
    finite_indices::Vector{Int},
    config::PreprocessingConfig,
    warnings::Vector{String},
)
    config.segment_policy == :all_finite && return finite_indices, nothing
    config.segment_policy == :recorded_ramp || return _preprocessing_error(
        experiment, "unsupported segment policy $(config.segment_policy)"
    )

    target_segment = experiment.kind == :dynamic ? 2 : 3
    if isnothing(experiment.segment)
        push!(
            warnings,
            "recorded segment is unavailable; used all finite rows before temperature cropping",
        )
        return finite_indices, nothing
    end
    selected = filter(finite_indices) do index
        value = experiment.segment[index]
        return !ismissing(value) && value == target_segment
    end
    isempty(selected) && _preprocessing_error(
        experiment, "recorded ramp segment $target_segment contains no finite rows"
    )
    return selected, target_segment
end

function _global_linear_heating_rate(
    experiment::Experiment, time_min::Vector{Float64}, temperature_K::Vector{Float64}
)
    centered_time = time_min .- mean(time_min)
    denominator = dot(centered_time, centered_time)
    denominator > 0 || _preprocessing_error(experiment, "time span is zero")
    rate = dot(centered_time, temperature_K .- mean(temperature_K)) / denominator
    isfinite(rate) && rate > 0 ||
        _preprocessing_error(experiment, "fitted heating rate is not finite and positive")
    return rate
end

"""
    cumulative_trapezoid(x, y; initial=0.0)

Cumulative trapezoidal integral along the supplied acquisition order. `x` need not be
uniform, but both vectors must be finite and have equal lengths.
"""
function cumulative_trapezoid(x::AbstractVector, y::AbstractVector; initial::Real=0.0)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal lengths"))
    isempty(x) && return Float64[]
    x_values = Float64.(x)
    y_values = Float64.(y)
    all(isfinite, x_values) && all(isfinite, y_values) ||
        throw(ArgumentError("x and y must contain only finite values"))
    integral = zeros(Float64, length(x_values))
    integral[1] = Float64(initial)
    for index in 2:length(x_values)
        integral[index] =
            integral[index - 1] +
            (x_values[index] - x_values[index - 1]) *
            (y_values[index] + y_values[index - 1]) / 2
    end
    return integral
end

"""
    preprocessing_fingerprint(config, initial_temperature_K, final_temperature_K)

Return a deterministic SHA-256 identifier for every numerical preprocessing choice and the
two mass-reference temperatures.
"""
function preprocessing_fingerprint(
    config::PreprocessingConfig, initial_temperature_K::Real, final_temperature_K::Real
)
    fields = (
        config.profile,
        config.invalid_row_policy,
        config.segment_policy,
        config.rebase_time,
        config.reference_mass_method,
        config.reference_half_window_K,
        config.smoothing_method,
        config.smoothing_half_window_K,
        config.local_polynomial_degree,
        config.derivative_method,
        config.heating_rate_method,
        config.monotonic_conversion_policy,
        config.analysis_conversion_range,
        config.reconstruction_rmse_tolerance,
        Float64(initial_temperature_K),
        Float64(final_temperature_K),
    )
    return bytes2hex(sha256(join(repr.(fields), "|")))
end

"""
    preprocess(experiment, analysis_config; initial_temperature_celsius, final_temperature_celsius,
               preprocessing_config)

Transform one raw experiment without mutating it. Invalid-row and segment selection precede
exact-temperature mass-reference estimation; interval endpoints are inclusive. Conversion is
never clipped. By default, local cubic polynomials in time use a physical temperature window
to smooth conversion and calculate both `dalpha/dt` and the local heating rate; their ratio is
`dalpha/dT`. Integration diagnostics always reconstruct the configured analysis curve from
`dalpha/dT`.
"""
function preprocess(
    experiment::Experiment,
    analysis_config::AnalysisConfig;
    initial_temperature_celsius::Real=analysis_config.conversion.initial_temperature_celsius,
    final_temperature_celsius::Real=analysis_config.conversion.final_temperature_celsius,
    preprocessing_config::PreprocessingConfig=analysis_config.preprocessing,
)
    initial_K = Float64(initial_temperature_celsius) + 273.15
    final_K = Float64(final_temperature_celsius) + 273.15
    initial_K < final_K || _preprocessing_error(
        experiment, "initial temperature must be below final temperature"
    )
    config = preprocessing_config
    warnings = String[]

    finite_mask = valid_row_mask(experiment)
    finite_indices = findall(finite_mask)
    invalid_count = length(finite_mask) - length(finite_indices)
    if invalid_count > 0 && config.invalid_row_policy == :error
        _preprocessing_error(
            experiment, "$invalid_count non-finite core row(s) are present"
        )
    elseif config.invalid_row_policy != :drop_nonfinite
        _preprocessing_error(
            experiment, "unsupported invalid-row policy $(config.invalid_row_policy)"
        )
    elseif invalid_count > 0
        push!(warnings, "dropped $invalid_count non-finite core row(s)")
    end

    segment_indices, selected_segment = _select_segment_rows(
        experiment, finite_indices, config, warnings
    )
    segment_temperature = experiment.temperature_K[segment_indices]
    segment_mass = experiment.mass_percent[segment_indices]
    minimum_temperature, maximum_temperature = extrema(segment_temperature)
    minimum_temperature <= initial_K < final_K <= maximum_temperature ||
        _preprocessing_error(
            experiment,
            "reference interval $(initial_temperature_celsius)–$(final_temperature_celsius) °C is not covered by the selected ramp ($(round(minimum_temperature - 273.15; digits=3))–$(round(maximum_temperature - 273.15; digits=3)) °C)",
        )

    initial_mass, initial_reference_count = _reference_mass(
        experiment, segment_temperature, segment_mass, initial_K, config
    )
    final_mass, final_reference_count = _reference_mass(
        experiment, segment_temperature, segment_mass, final_K, config
    )
    mass_loss = initial_mass - final_mass
    isfinite(mass_loss) && mass_loss > 0 || _preprocessing_error(
        experiment,
        "initial reference mass ($initial_mass %) must exceed final reference mass ($final_mass %)",
    )

    retained_indices = filter(segment_indices) do index
        return initial_K <= experiment.temperature_K[index] <= final_K
    end
    minimum_points = max(3, config.local_polynomial_degree + 2)
    length(retained_indices) >= minimum_points || _preprocessing_error(
        experiment, "analysis interval contains fewer than $minimum_points rows"
    )

    temperature = copy(experiment.temperature_K[retained_indices])
    raw_time = copy(experiment.time_min[retained_indices])
    all(>(0), diff(raw_time)) || _preprocessing_error(
        experiment, "retained acquisition time is not strictly increasing"
    )
    time = config.rebase_time ? raw_time .- first(raw_time) : copy(raw_time)
    mass = copy(experiment.mass_percent[retained_indices])
    mass_fraction = mass ./ initial_mass
    alpha = (initial_mass .- mass) ./ mass_loss

    needs_local_fit =
        config.smoothing_method == :local_polynomial ||
        config.derivative_method == :local_polynomial
    local_fit, local_derivative = if needs_local_fit
        local_polynomial_estimate(
            time,
            alpha,
            temperature;
            half_window=config.smoothing_half_window_K,
            degree=config.local_polynomial_degree,
        )
    else
        copy(alpha), Float64[]
    end

    analysis_alpha = if config.smoothing_method == :none
        copy(alpha)
    elseif config.smoothing_method == :local_polynomial
        local_fit
    else
        _preprocessing_error(
            experiment, "unsupported smoothing method $(config.smoothing_method)"
        )
    end
    dalpha_dt = if config.derivative_method == :finite_difference
        finite_difference_derivative(time, analysis_alpha)
    elseif config.derivative_method == :local_polynomial
        local_derivative
    else
        _preprocessing_error(
            experiment, "unsupported derivative method $(config.derivative_method)"
        )
    end

    fitted_rate = _global_linear_heating_rate(experiment, time, temperature)
    heating_rate = if config.heating_rate_method == :global_linear
        fill(fitted_rate, length(time))
    elseif config.heating_rate_method == :local_polynomial
        _, local_rate = local_polynomial_estimate(
            time,
            temperature,
            temperature;
            half_window=config.smoothing_half_window_K,
            degree=config.local_polynomial_degree,
        )
        all(value -> isfinite(value) && value > 0, local_rate) || _preprocessing_error(
            experiment,
            "local-polynomial heating-rate curve is not finite and positive",
        )
        local_rate
    else
        _preprocessing_error(
            experiment, "unsupported heating-rate method $(config.heating_rate_method)"
        )
    end
    dalpha_dT = dalpha_dt ./ heating_rate
    reconstructed = cumulative_trapezoid(temperature, dalpha_dT; initial=analysis_alpha[1])
    residual = reconstructed .- analysis_alpha
    reconstruction_rmse = sqrt(mean(abs2, residual))
    reconstruction_max_error = maximum(abs, residual)

    unit_interval_tolerance = 1.0e-8
    outside_count = count(
        value -> value < -unit_interval_tolerance || value > 1 + unit_interval_tolerance,
        alpha,
    )
    reversal_count = count(<(-1.0e-8), diff(analysis_alpha))
    negative_derivative_count = count(<(-1.0e-10), dalpha_dT)
    nonincreasing_temperature_count = count(<=(0), diff(temperature))
    outside_count > 0 && push!(
        warnings,
        "$outside_count conversion value(s) lie outside [0, 1]; values were not clipped",
    )
    reversal_count > 0 && push!(
        warnings,
        "$reversal_count analysis-conversion step(s) decrease by more than 1e-8",
    )
    negative_derivative_count > 0 && push!(
        warnings, "$negative_derivative_count derivative value(s) are below -1e-10 K^-1"
    )
    nonincreasing_temperature_count > 0 && push!(
        warnings,
        "$nonincreasing_temperature_count retained temperature step(s) are zero or negative",
    )
    reconstruction_rmse > config.reconstruction_rmse_tolerance && push!(
        warnings,
        "derivative reconstruction RMSE $(round(reconstruction_rmse; sigdigits=5)) exceeds $(config.reconstruction_rmse_tolerance)",
    )
    if config.monotonic_conversion_policy == :error && reversal_count > 0
        _preprocessing_error(
            experiment,
            "$reversal_count conversion reversals violate monotonic_conversion_policy=error",
        )
    elseif config.monotonic_conversion_policy != :diagnose &&
        config.monotonic_conversion_policy != :error
        _preprocessing_error(
            experiment,
            "unsupported monotonic-conversion policy $(config.monotonic_conversion_policy)",
        )
    end

    diagnostics = PreprocessingDiagnostics(
        length(finite_mask),
        length(finite_indices),
        length(segment_indices),
        length(retained_indices),
        invalid_count,
        length(finite_indices) - length(segment_indices),
        length(segment_indices) - length(retained_indices),
        selected_segment,
        initial_reference_count,
        final_reference_count,
        initial_mass,
        final_mass,
        fitted_rate,
        minimum(alpha),
        maximum(alpha),
        outside_count,
        reversal_count,
        negative_derivative_count,
        reconstruction_rmse,
        reconstruction_max_error,
        warnings,
    )

    return ProcessedExperiment(
        experiment.id,
        experiment.source_file,
        experiment.source_variable,
        experiment.source_sha256,
        experiment.composition,
        copy(retained_indices),
        temperature,
        raw_time,
        time,
        mass,
        mass_fraction,
        alpha,
        analysis_alpha,
        dalpha_dT,
        dalpha_dt,
        heating_rate,
        reconstructed,
        residual,
        initial_K,
        final_K,
        config,
        preprocessing_fingerprint(config, initial_K, final_K),
        diagnostics,
    )
end
