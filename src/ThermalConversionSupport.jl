"""
An acquisition-order encounter with a target conversion. `kind` is `:upward`,
`:downward`, `:tangency`, `:endpoint`, or `:plateau`. Segment-local rows and
physical source-line bounds identify the recorded bracket. A plateau reports its
full measured temperature range and has no single crossing temperature or time.
An isolated exact sample uses a one-row bracket. Time is reported in the original
clock origin, in minutes, when recorded.
"""
struct TGConversionCrossing
    kind::Symbol
    segment_rows::UnitRange{Int}
    source_lines::Tuple{Int,Int}
    temperature_K::Union{Nothing,Float64}
    time_min::Union{Nothing,Float64}
    temperature_range_K::Tuple{Float64,Float64}
end

"""
All encounters with `alpha` in one TG curve, in acquisition order. `status` is
`:missing` for no encounter and `:ambiguous` for multiple encounters. A single
encounter has status `:unique_upward`, `:unique_downward`, `:tangency`,
`:endpoint`, or `:plateau`.
"""
struct TGCrossingDiagnostic
    alpha::Float64
    status::Symbol
    crossings::Vector{TGConversionCrossing}
end

"""
Group coverage at one target conversion. `crossings` follows the group run order.
`common_upward` requires a unique interior upward crossing in every run.
`upward_run_count` counts the individual runs that meet that criterion. At common
points, temperature statistics summarize the interpolated crossing temperatures.
The reciprocal span and centered sum of squares use `1/T` in inverse Kelvin.
All temperature summaries are `nothing` when a common crossing is unavailable.
"""
struct TGConversionSupportPoint
    alpha::Float64
    crossings::Vector{TGCrossingDiagnostic}
    common_upward::Bool
    upward_run_count::Int
    temperature_min_K::Union{Nothing,Float64}
    temperature_max_K::Union{Nothing,Float64}
    temperature_mean_K::Union{Nothing,Float64}
    temperature_span_K::Union{Nothing,Float64}
    reciprocal_temperature_span_K_inv::Union{Nothing,Float64}
    reciprocal_temperature_sxx_K_inv2::Union{Nothing,Float64}
end

"""
Conversion support for a declared set of TG runs. `run_ids` preserves the input
order and `missing_time_run_ids` identifies runs without a recorded clock.
`mass_basis` is `:absolute_mass` or `:relative_mass`. Each entry in `points`
describes coverage at one requested conversion.
"""
struct TGConversionGroupSupport
    material_id::String
    mass_basis::Symbol
    run_ids::Vector{String}
    missing_time_run_ids::Vector{String}
    points::Vector{TGConversionSupportPoint}
end

function _validate_tg_support_curve(curve)
    n = length(curve.alpha)
    n > 0 &&
    n ==
    length(curve.temperature_K) ==
    length(curve.source_lines) ==
    length(curve.reaction_rows) ||
        throw(ArgumentError("$(curve.run_id) has inconsistent TG support rows"))
    all(isfinite, curve.alpha) ||
        throw(ArgumentError("$(curve.run_id) has nonfinite conversion"))
    all(T -> isfinite(T) && T > 0, curve.temperature_K) ||
        throw(ArgumentError("$(curve.run_id) needs finite positive temperatures in Kelvin"))
    first(curve.reaction_rows) >= 1 &&
    first(curve.source_lines) >= 1 &&
    all(>(0), diff(curve.source_lines)) ||
        throw(ArgumentError("$(curve.run_id) needs ordered physical source lines"))
    t = curve.time_min
    if t !== nothing
        length(t) == n && all(isfinite, t) && all(>(0), diff(t)) ||
            throw(ArgumentError("$(curve.run_id) needs finite strictly increasing time"))
    end
    return nothing
end

function _tg_support_target(alpha)
    isfinite(alpha) && 0 <= alpha <= 1 ||
        throw(ArgumentError("target conversion must be finite and within [0, 1]"))
    return Float64(alpha)
end

function _tg_level_event(curve, kind, first_index, last_index, fraction)
    temperature = curve.temperature_K
    clock = curve.time_min
    temp_range = extrema(view(temperature, first_index:last_index))
    at(values) = (1 - fraction) * values[first_index] + fraction * values[last_index]
    return TGConversionCrossing(
        kind,
        curve.reaction_rows[first_index:last_index],
        (curve.source_lines[first_index], curve.source_lines[last_index]),
        fraction === nothing ? nothing : at(temperature),
        fraction === nothing || clock === nothing ? nothing : at(clock),
        temp_range,
    )
end

function _tg_exact_level_kind(alpha, target, i, j)
    i != j && return :plateau
    (i == firstindex(alpha) || i == lastindex(alpha)) && return :endpoint
    below_before = alpha[i - 1] < target
    below_after = alpha[i + 1] < target
    below_before == below_after && return :tangency
    return below_before ? :upward : :downward
end

function _tg_crossing_fraction(left, right, target)
    width = right - left
    # Finite observations of opposite signs can still have an overflowing difference.
    if isfinite(width)
        return (target - left) / width
    end
    return (target / 2 - left / 2) / (right / 2 - left / 2)
end

function _tg_conversion_crossings(curve, target)
    alpha = curve.alpha
    events = TGConversionCrossing[]
    i = firstindex(alpha)
    while i <= lastindex(alpha)
        if alpha[i] == target
            j = i
            while j < lastindex(alpha) && alpha[j + 1] == target
                j += 1
            end
            kind = _tg_exact_level_kind(alpha, target, i, j)
            fraction = kind === :plateau ? nothing : 0.0
            push!(events, _tg_level_event(curve, kind, i, j, fraction))
            i = j + 1
            continue
        end
        if i < lastindex(alpha) &&
            alpha[i + 1] != target &&
            (alpha[i] < target) != (alpha[i + 1] < target)
            kind = alpha[i] < target ? :upward : :downward
            fraction = _tg_crossing_fraction(alpha[i], alpha[i + 1], target)
            push!(events, _tg_level_event(curve, kind, i, i + 1, fraction))
        end
        i += 1
    end
    status = if isempty(events)
        :missing
    elseif length(events) > 1
        :ambiguous
    else
        kind = only(events).kind
        if kind === :upward
            :unique_upward
        elseif kind === :downward
            :unique_downward
        else
            kind
        end
    end
    return TGCrossingDiagnostic(target, status, events)
end

"""
    tg_conversion_crossings(curve, alpha)

Find all encounters with a target conversion in a [`TGConversionResult`](@ref).
Strict brackets use acquisition-order linear interpolation of measured temperature
and optional time. Exact target comparisons use the recorded numerical values.
Consecutive exact samples are one plateau. A single exact interior sample is a
crossing when its neighbors straddle the target and a tangency when they lie on
the same side. An exact boundary sample is a one-sided endpoint encounter.
"""
function tg_conversion_crossings(curve::TGConversionResult, alpha::Real)
    target = _tg_support_target(alpha)
    _validate_tg_support_curve(curve)
    return _tg_conversion_crossings(curve, target)
end

function _tg_support_material(curves, material_id)
    if material_id === nothing
        material = first(curves).material_id
        material === nothing &&
            throw(ArgumentError("declare material_id for TG runs with unknown material"))
    else
        material = String(strip(material_id))
        isempty(material) && throw(ArgumentError("material_id cannot be empty"))
    end
    for curve in curves
        known = curve.material_id
        if known === nothing
            material_id === nothing && throw(
                ArgumentError("declare material_id for TG runs with unknown material")
            )
        elseif known != material
            throw(ArgumentError("TG run $(curve.run_id) belongs to a different material"))
        end
    end
    return material
end

function _tg_support_point(curves, target)
    diagnostics = [_tg_conversion_crossings(curve, target) for curve in curves]
    upward_count = count(d -> d.status === :unique_upward, diagnostics)
    if upward_count != length(curves)
        return TGConversionSupportPoint(
            target,
            diagnostics,
            false,
            upward_count,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
        )
    end
    temperatures = Float64[only(d.crossings).temperature_K for d in diagnostics]
    lower, upper = extrema(temperatures)
    mean_temperature = mean(temperatures)
    inverse_T = inv.(temperatures)
    all(isfinite, inverse_T) && isfinite(mean_temperature) ||
        throw(ArgumentError("temperature statistics overflow at conversion $target"))
    inv_lower, inv_upper = extrema(inverse_T)
    inv_center = mean(inverse_T)
    sxx = sum(x -> abs2(x - inv_center), inverse_T)
    isfinite(sxx) || throw(
        ArgumentError(
            "reciprocal temperature sum of squares overflows at conversion $target"
        ),
    )
    return TGConversionSupportPoint(
        target,
        diagnostics,
        true,
        upward_count,
        lower,
        upper,
        mean_temperature,
        upper - lower,
        inv_upper - inv_lower,
        sxx,
    )
end

"""
    tg_conversion_support(curves, targets; material_id=nothing)

Check a strictly increasing target grid on a nonempty group of TG runs. Runs need
distinct identifiers and source fingerprints and must describe the same material.
Supply `material_id` to declare that runs with unknown labels belong to that material.
Known labels must agree with the group. All mass channels must use the same absolute
or relative mass basis. Common coverage requires a unique interior upward crossing
per run. A recorded clock is optional for temperature coverage.

Temperature and reciprocal-temperature spans describe coverage. Assess the
conversion definition, temperature programs, and suitability of the coverage for
the intended kinetic method separately.
"""
function tg_conversion_support(
    curves::AbstractVector{TGConversionResult},
    targets::AbstractVector{<:Real};
    material_id::Union{Nothing,AbstractString}=nothing,
)
    isempty(curves) && throw(ArgumentError("TG support needs at least one run"))
    isempty(targets) && throw(ArgumentError("TG support needs at least one target"))
    grid = _tg_support_target.(targets)
    all(>(0), diff(grid)) ||
        throw(ArgumentError("TG support targets must be strictly increasing"))
    material = _tg_support_material(curves, material_id)
    run_ids = [curve.run_id for curve in curves]
    allunique(run_ids) || throw(ArgumentError("TG support needs distinct run IDs"))
    allunique(curve.source_sha256 for curve in curves) ||
        throw(ArgumentError("TG support needs distinct source fingerprints"))
    basis(unit) = first(get(_THERMAL_SIGNAL_UNITS, unit, (:unmapped, 1.0)))
    mass_basis = basis(first(curves).mass_unit)
    mass_basis in (:absolute_mass, :relative_mass) &&
    all(curve -> basis(curve.mass_unit) === mass_basis, curves) ||
        throw(ArgumentError("TG runs must use the same absolute or relative mass basis"))
    foreach(_validate_tg_support_curve, curves)
    missing_time = [curve.run_id for curve in curves if curve.time_min === nothing]
    points = [_tg_support_point(curves, target) for target in grid]
    return TGConversionGroupSupport(material, mass_basis, run_ids, missing_time, points)
end
