"""
One requested TG conversion increment and its recorded-time branch diagnostic.
`status === :supported` requires unique interior upward boundary crossings,
recorded time, and strictly increasing conversion throughout the clipped interval.
The measured interior nodes are unchanged. A boundary node has the source-line
bracket used to interpolate it, or equal source lines for an exact sample.

Path arrays are empty when the boundaries or clock are unavailable. Reversed or
flat paths retain their nodes, backtracking, and plateau duration for inspection.
`max_backtracking_alpha` is the largest drop below any earlier conversion in the
clipped interval. The original clock is in minutes.
"""
struct TGRecordedTimeInterval
    run_id::String
    segment_index::Int
    source_sha256::String
    alpha_lower::Float64
    alpha_upper::Float64
    lower_crossing_status::Symbol
    upper_crossing_status::Symbol
    status::Symbol
    alpha::Vector{Float64}
    time_min::Vector{Float64}
    temperature_K::Vector{Float64}
    source_lines::Vector{Tuple{Int,Int}}
    reversal_steps::Union{Nothing,Int}
    plateau_steps::Union{Nothing,Int}
    max_backtracking_alpha::Union{Nothing,Float64}
    plateau_duration_min::Union{Nothing,Float64}
end

"""
Logarithms of four Arrhenius integrals at one declared probe energy. All integrals
use minutes. `recorded_log_integral` uses the clipped measured nodes.
`legacy_log_integral` uses uniform-conversion sampling on the same branch.
`refined_log_integral` and `doubled_log_integral` use equal time subdivisions of
each measured segment and piecewise-linear measured temperature.
"""
struct TGRecordedTimeQuadrature
    probe_energy_kJ_per_mol::Float64
    legacy_alpha_nodes::Int
    subdivisions_per_step::Int
    recorded_log_integral::Float64
    legacy_log_integral::Float64
    refined_log_integral::Float64
    doubled_log_integral::Float64
end

function _tg_unavailable_time_interval(
    curve, lower, upper, lower_status, upper_status, status
)
    return TGRecordedTimeInterval(
        curve.run_id,
        curve.segment_index,
        curve.source_sha256,
        lower,
        upper,
        lower_status,
        upper_status,
        status,
        Float64[],
        Float64[],
        Float64[],
        Tuple{Int,Int}[],
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function _tg_time_path(curve, lower, upper, start, stop)
    offset = first(curve.reaction_rows) - 1
    # Exact boundary rows are supplied by the crossing, so keep them out of the interior.
    inner = (first(start.segment_rows) - offset + 1):(last(stop.segment_rows) - offset - 1)
    alpha = [lower; curve.alpha[inner]; upper]
    t = [start.time_min::Float64; curve.time_min[inner]; stop.time_min::Float64]
    T = [
        start.temperature_K::Float64;
        curve.temperature_K[inner];
        stop.temperature_K::Float64
    ]
    lines = [(line, line) for line in curve.source_lines[inner]]
    pushfirst!(lines, start.source_lines)
    push!(lines, stop.source_lines)
    return alpha, t, T, lines
end

"""
    tg_recorded_time_interval(curve, alpha; delta_alpha)

Check the complete recorded branch over `[alpha - delta_alpha, alpha]` on one
[`TGConversionResult`](@ref). Both targets must be strictly inside (0, 1) and have
one interior upward crossing in the selected conversion record. Only the boundary
temperature and time are interpolated. Every recorded interior row remains in
acquisition order.

An interior decrease gives `:branch_reversal`. An exact flat step gives
`:branch_plateau` when there are no reversals. Both block strict support. The path
diagnostics preserve the size of these features without applying a noise tolerance.
Other unavailable statuses describe the target bounds, boundary crossings, or time.
"""
function tg_recorded_time_interval(
    curve::TGConversionResult, alpha::Real; delta_alpha::Real
)
    upper = _tg_support_target(alpha)
    isfinite(delta_alpha) && 0 < delta_alpha < 1 ||
        throw(ArgumentError("delta_alpha must lie strictly between zero and one"))
    lower = upper - Float64(delta_alpha)
    if !(0 < lower < upper < 1)
        return _tg_unavailable_time_interval(
            curve, lower, upper, :unchecked, :unchecked, :outside_interior
        )
    end
    _validate_tg_support_curve(curve)
    start = _tg_conversion_crossings(curve, lower)
    stop = _tg_conversion_crossings(curve, upper)
    if start.status !== :unique_upward || stop.status !== :unique_upward
        return _tg_unavailable_time_interval(
            curve, lower, upper, start.status, stop.status, :boundary_not_unique
        )
    end
    if curve.time_min === nothing
        return _tg_unavailable_time_interval(
            curve, lower, upper, start.status, stop.status, :missing_time
        )
    end

    path_alpha, t, T, lines = _tg_time_path(
        curve, lower, upper, only(start.crossings), only(stop.crossings)
    )
    dt = diff(t)
    if !all(step -> isfinite(step) && step > 0, dt) || !isfinite(last(t) - first(t))
        return _tg_unavailable_time_interval(
            curve, lower, upper, start.status, stop.status, :clock_not_increasing
        )
    end
    steps = diff(path_alpha)
    reversals = count(<(0), steps)
    plateaus = count(iszero, steps)
    backtrack = maximum(accumulate(max, path_alpha) .- path_alpha)
    flat_duration = sum(dt[iszero.(steps)])
    status = if reversals > 0
        :branch_reversal
    elseif plateaus > 0
        :branch_plateau
    else
        :supported
    end
    return TGRecordedTimeInterval(
        curve.run_id,
        curve.segment_index,
        curve.source_sha256,
        lower,
        upper,
        start.status,
        stop.status,
        status,
        path_alpha,
        t,
        T,
        lines,
        reversals,
        plateaus,
        backtrack,
        flat_duration,
    )
end

function _tg_uniform_alpha_time_path(interval, nodes)
    alpha, t, T = interval.alpha, interval.time_min, interval.temperature_K
    resampled_t = Vector{Float64}(undef, nodes)
    resampled_T = similar(resampled_t)
    for (j, target) in enumerate(range(first(alpha), last(alpha); length=nodes))
        i = min(searchsortedlast(alpha, target), length(alpha) - 1)
        fraction = (target - alpha[i]) / (alpha[i + 1] - alpha[i])
        resampled_t[j] = (1 - fraction) * t[i] + fraction * t[i + 1]
        resampled_T[j] = (1 - fraction) * T[i] + fraction * T[i + 1]
    end
    return resampled_t, resampled_T
end

function _tg_log_time_subdivided_integral(E, t, T, subdivisions)
    reference_temperature = E >= 0 ? maximum(T) : minimum(T)
    shift = -E / (GAS_CONSTANT_KJ_PER_MOL_K * reference_temperature)
    scaled_rate(temp) = exp(-E / (GAS_CONSTANT_KJ_PER_MOL_K * temp) - shift)
    total = 0.0
    left = scaled_rate(first(T))
    for i in 2:length(t)
        right = scaled_rate(T[i])
        segment = (left + right) / 2
        for k in 1:(subdivisions - 1)
            fraction = k / subdivisions
            temp = (1 - fraction) * T[i - 1] + fraction * T[i]
            segment += scaled_rate(temp)
        end
        total += (t[i] - t[i - 1]) * (segment / subdivisions)
        left = right
    end
    return shift + log(total)
end

"""
    tg_recorded_time_quadrature(interval, probe_energy_kJ_per_mol;
                                legacy_alpha_nodes=41, subdivisions_per_step=8)

Compare trapezoid integrals over a strictly supported
[`TGRecordedTimeInterval`](@ref) at a fixed, nonnegative probe energy.
The legacy comparison uses uniform conversion sampling on the same branch.
The refined comparisons subdivide each time step while keeping the measured
piecewise-linear temperature history. Their difference measures numerical
refinement on that history. It cannot resolve the temperature between observations.
"""
function tg_recorded_time_quadrature(
    interval::TGRecordedTimeInterval,
    probe_energy_kJ_per_mol::Real;
    legacy_alpha_nodes::Integer=41,
    subdivisions_per_step::Integer=8,
)
    interval.status === :supported || throw(
        ArgumentError("recorded-time quadrature requires a supported conversion interval"),
    )
    energy = Float64(probe_energy_kJ_per_mol)
    isfinite(energy) && energy >= 0 ||
        throw(ArgumentError("probe energy must be finite and nonnegative"))
    legacy_alpha_nodes >= 3 ||
        throw(ArgumentError("legacy conversion sampling needs at least three nodes"))
    2 <= subdivisions_per_step <= typemax(Int) ÷ 2 ||
        throw(ArgumentError("time refinement requires at least two subdivisions per step"))
    t, T = interval.time_min, interval.temperature_K
    legacy_t, legacy_T = _tg_uniform_alpha_time_path(interval, legacy_alpha_nodes)
    all(>(0), diff(legacy_t)) ||
        throw(ArgumentError("legacy conversion sampling exceeds recorded clock precision"))

    native = _tg_log_time_subdivided_integral(energy, t, T, 1)
    legacy = _tg_log_time_subdivided_integral(energy, legacy_t, legacy_T, 1)
    refined = _tg_log_time_subdivided_integral(energy, t, T, subdivisions_per_step)
    doubled = _tg_log_time_subdivided_integral(energy, t, T, 2 * subdivisions_per_step)
    all(isfinite, (native, legacy, refined, doubled)) ||
        throw(ArgumentError("Arrhenius integral is outside numerical range"))
    return TGRecordedTimeQuadrature(
        energy,
        Int(legacy_alpha_nodes),
        Int(subdivisions_per_step),
        native,
        legacy,
        refined,
        doubled,
    )
end
