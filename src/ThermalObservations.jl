"""
    ThermalColumnSpec(name, unit)

An explicit source column name and unit. No units are inferred from headers.
Temperature mappings accept `:K` and `:C`; time mappings accept `:s`, `:min`, and `:h`.
"""
struct ThermalColumnSpec
    name::String
    unit::Symbol

    function ThermalColumnSpec(name::AbstractString, unit::Symbol)
        label = String(strip(name))
        isempty(label) && throw(ArgumentError("column name cannot be empty"))
        isempty(String(unit)) && throw(ArgumentError("column unit cannot be empty"))
        return new(label, unit)
    end
end

const _THERMAL_SIGNAL_UNITS = Dict(
    :kg => (:absolute_mass, 1.0),
    :g => (:absolute_mass, 1e-3),
    :mg => (:absolute_mass, 1e-6),
    :ug => (:absolute_mass, 1e-9),
    :fraction => (:relative_mass, 1.0),
    :percent => (:relative_mass, 1e-2),
    :W => (:heat_flow, 1.0),
    :mW => (:heat_flow, 1e-3),
    :uW => (:heat_flow, 1e-6),
    :W_per_kg => (:specific_heat_flow, 1.0),
    :W_per_g => (:specific_heat_flow, 1e3),
    :mW_per_mg => (:specific_heat_flow, 1e3),
    :mW_per_g => (:specific_heat_flow, 1.0),
)

"""
    ThermalSignalSpec(id, column, quantity; exotherm_sign=nothing)

A measured channel's identifier, source column, and physical quantity. Standard quantities
include `:mass`, `:mass_loss`, `:heat_flow`, and `:specific_heat_flow`. Other tagged
quantities can be retained in their source units. For heat flow, `exotherm_sign` is `-1`
or `1` when known. Recording it never flips values.
"""
struct ThermalSignalSpec
    id::Symbol
    column::ThermalColumnSpec
    quantity::Symbol
    exotherm_sign::Union{Nothing,Int}

    function ThermalSignalSpec(
        id::Symbol, column::ThermalColumnSpec, quantity::Symbol; exotherm_sign=nothing
    )
        isempty(String(id)) && throw(ArgumentError("signal identifier cannot be empty"))
        isempty(String(quantity)) && throw(ArgumentError("signal quantity cannot be empty"))
        dimension = first(get(_THERMAL_SIGNAL_UNITS, column.unit, (:unmapped, 1.0)))
        if quantity in (:mass, :mass_loss)
            dimension in (:absolute_mass, :relative_mass) ||
                throw(ArgumentError("invalid mass unit $(column.unit) for $id"))
        elseif quantity in (:heat_flow, :specific_heat_flow)
            dimension == quantity ||
                throw(ArgumentError("unit $(column.unit) is incompatible with $quantity"))
        end
        sign = if exotherm_sign === nothing
            nothing
        else
            quantity in (:heat_flow, :specific_heat_flow) ||
                throw(ArgumentError("exotherm sign applies only to heat flow"))
            exotherm_sign isa Integer &&
            !(exotherm_sign isa Bool) &&
            exotherm_sign in (-1, 1) ||
                throw(ArgumentError("exotherm sign must be -1, 1, or nothing"))
            Int(exotherm_sign)
        end
        return new(id, column, quantity, sign)
    end
end

"""
    ThermalColumnMapping(temperature, signals; time=nothing, segment_column=nothing)

Select distinct source columns with declared units. Use `time=nothing` for sources
without recorded time and `segment_column=nothing` for a single unlabelled block.
Import does not infer time, heating rate, or physical program phases.
"""
struct ThermalColumnMapping
    temperature::ThermalColumnSpec
    signals::Vector{ThermalSignalSpec}
    time::Union{Nothing,ThermalColumnSpec}
    segment_column::Union{Nothing,String}

    function ThermalColumnMapping(
        temperature::ThermalColumnSpec, signals; time=nothing, segment_column=nothing
    )
        temperature.unit in (:K, :C) ||
            throw(ArgumentError("temperature unit must be :K or :C"))
        if time !== nothing
            time isa ThermalColumnSpec && time.unit in (:s, :min, :h) ||
                throw(ArgumentError("time must map a column in :s, :min, or :h"))
        end
        columns = ThermalSignalSpec[s for s in signals]
        isempty(columns) && throw(ArgumentError("map at least one signal"))
        allunique(s.id for s in columns) ||
            throw(ArgumentError("duplicate signal identifier"))
        segment = if segment_column === nothing
            nothing
        else
            label = String(strip(segment_column))
            isempty(label) && throw(ArgumentError("segment column cannot be empty"))
            label
        end
        names = [temperature.name; [s.column.name for s in columns]]
        time === nothing || push!(names, time.name)
        segment === nothing || push!(names, segment)
        allunique(names) || throw(ArgumentError("a source column cannot be mapped twice"))
        return new(temperature, columns, time, segment)
    end
end

"""
Measured values of one explicitly mapped channel, stored in acquisition order and
in the original unit. See [`signal_values`](@ref) for checked unit conversion.
"""
struct ThermalSignal
    spec::ThermalSignalSpec
    values::Vector{Float64}
end

"""
A contiguous acquisition block. `index` is its one-based ordinal in the original
run; `label` is its recorded string (or `nothing`). The same label can recur in
separate blocks. `rows` indexes the containing run, not physical source lines.
"""
struct ThermalSegment
    index::Int
    label::Union{Nothing,String}
    rows::UnitRange{Int}
end

"""
Import source identity: the original path, byte SHA-256, physical source line for
each numerical row, and headers of unselected columns. A segment selection retains
the original hash and corresponding line numbers.
"""
struct ThermalProvenance
    source_path::String
    sha256::String
    source_lines::Vector{Int}
    unmapped_columns::Vector{String}
end

"""
A material-independent thermal run returned by a validated importer. Axis and channel
vectors hold source values, not normalized conversion. Temperature is mandatory;
recorded time and material identity are optional. Labels, segment boundaries, metadata,
and source provenance are retained.
"""
struct ThermalRun
    id::String
    material_id::Union{Nothing,String}
    temperature::Vector{Float64}
    temperature_unit::Symbol
    time::Union{Nothing,Vector{Float64}}
    time_unit::Union{Nothing,Symbol}
    signals::Vector{ThermalSignal}
    segments::Vector{ThermalSegment}
    metadata::Dict{String,String}
    provenance::ThermalProvenance
end

_thermal_temperature_K(x, unit) = unit === :C ? x + 273.15 : x
_thermal_time_factor(unit) =
    if unit === :s
        1 / 60
    elseif unit === :h
        60.0
    else
        1.0
    end

function _thermal_segments(n, labels)
    labels === nothing && return [ThermalSegment(1, nothing, 1:n)]
    segments = ThermalSegment[]
    first_row = 1
    for i in 2:n
        if labels[i] != labels[i - 1]
            push!(
                segments,
                ThermalSegment(length(segments) + 1, labels[i - 1], first_row:(i - 1)),
            )
            first_row = i
        end
    end
    push!(segments, ThermalSegment(length(segments) + 1, labels[end], first_row:n))
    return segments
end

function _build_thermal_run(
    id,
    material_id,
    mapping,
    raw_temperature,
    raw_time,
    raw_signals,
    labels,
    metadata,
    provenance,
)
    fail(message) = throw(DataImportError(provenance.source_path, "", message))
    run_id = String(strip(id))
    isempty(run_id) && throw(ArgumentError("thermal run id cannot be empty"))
    material = material_id === nothing ? nothing : String(strip(material_id))
    material == "" &&
        throw(ArgumentError("material id cannot be empty; use nothing when unknown"))
    n = length(raw_temperature)
    n > 0 || fail("no numeric data rows")
    length(provenance.source_lines) == n || fail("source line count does not match data")
    if mapping.time === nothing
        raw_time === nothing || fail("time values supplied without a time mapping")
    else
        raw_time !== nothing && length(raw_time) == n ||
            fail("time length does not match data")
    end
    length(raw_signals) == length(mapping.signals) ||
        fail("mapped signal count does not match data")
    all(length(values) == n for values in raw_signals) ||
        fail("signal length does not match data")
    if mapping.segment_column === nothing
        labels === nothing || fail("segment labels supplied without a segment mapping")
    else
        labels !== nothing && length(labels) == n ||
            fail("segment labels do not match data")
        bad_label = findfirst(isempty, labels)
        bad_label === nothing ||
            fail("empty segment label at source line $(provenance.source_lines[bad_label])")
    end
    for i in eachindex(raw_temperature)
        kelvin = _thermal_temperature_K(raw_temperature[i], mapping.temperature.unit)
        isfinite(kelvin) && kelvin > 0 || fail(
            "temperature must be finite and above absolute zero at source line $(provenance.source_lines[i])",
        )
    end
    if raw_time !== nothing
        factor = _thermal_time_factor(mapping.time.unit)
        bad_time = findfirst(t -> !isfinite(t * factor), raw_time)
        bad_time === nothing ||
            fail("time must be finite at source line $(provenance.source_lines[bad_time])")
    end
    for (spec, values) in zip(mapping.signals, raw_signals)
        bad_signal = findfirst(x -> !isfinite(x), values)
        bad_signal === nothing || fail(
            "signal $(spec.id) must be finite at source line $(provenance.source_lines[bad_signal])",
        )
    end
    segments = _thermal_segments(n, labels)
    if raw_time !== nothing
        for segment in segments
            for i in (first(segment.rows) + 1):last(segment.rows)
                raw_time[i] > raw_time[i - 1] || fail(
                    "time must strictly increase within segment $(segment.index) at source line $(provenance.source_lines[i])",
                )
            end
        end
    end
    return ThermalRun(
        run_id,
        material,
        raw_temperature,
        mapping.temperature.unit,
        raw_time,
        mapping.time === nothing ? nothing : mapping.time.unit,
        [
            ThermalSignal(spec, values) for
            (spec, values) in zip(mapping.signals, raw_signals)
        ],
        segments,
        metadata,
        provenance,
    )
end

"""
    temperature_K(run)

Return a new vector of absolute temperatures in Kelvin, preserving acquisition order.
"""
function temperature_K(run::ThermalRun)
    return [_thermal_temperature_K(x, run.temperature_unit) for x in run.temperature]
end

"""
    time_min(run)

Return recorded times in minutes, or `nothing` if the source has no time axis.
The original clock origin and any between-segment resets are retained.
"""
function time_min(run::ThermalRun)
    return run.time === nothing ? nothing : run.time .* _thermal_time_factor(run.time_unit)
end

"""
    signal_values(run, id; unit=nothing)

Return a copy of a measured signal, optionally converted between compatible units.
Conversion never changes sign or uses a sample reference to convert between absolute
and relative mass, or between heat flow and specific heat flow. Unrecognized or
auxiliary quantity tags may be read in their source unit only.
"""
function signal_values(run::ThermalRun, id::Symbol; unit=nothing)
    index = findfirst(s -> s.spec.id === id, run.signals)
    index === nothing && throw(ArgumentError("run $(run.id) has no signal $id"))
    signal = run.signals[index]
    origin = signal.spec.column.unit
    (unit === nothing || unit === origin) && return copy(signal.values)
    signal.spec.quantity in (:mass, :mass_loss, :heat_flow, :specific_heat_flow) ||
        throw(ArgumentError("no conversion is defined for $(signal.spec.quantity)"))
    source_dimension, source_scale = _THERMAL_SIGNAL_UNITS[origin]
    destination_dimension, destination_scale = get(
        _THERMAL_SIGNAL_UNITS, unit, (:unmapped, 1.0)
    )
    source_dimension == destination_dimension || throw(
        ArgumentError(
            "cannot convert $id from $origin to $unit without another observation model"
        ),
    )
    values = signal.values .* (source_scale / destination_scale)
    all(isfinite, values) || throw(ArgumentError("signal unit conversion overflow for $id"))
    return values
end

"""
    select_thermal_segment(run, index)

Return an independent selection by original segment ordinal, not by label or a
guessed phase. Preserve its source label, recorded clock, signal units and line numbers.
"""
function select_thermal_segment(run::ThermalRun, index::Integer)
    i = findfirst(segment -> segment.index == index, run.segments)
    i === nothing && throw(ArgumentError("run $(run.id) has no segment ordinal $index"))
    block = run.segments[i]
    rows = block.rows
    provenance = ThermalProvenance(
        run.provenance.source_path,
        run.provenance.sha256,
        run.provenance.source_lines[rows],
        copy(run.provenance.unmapped_columns),
    )
    return ThermalRun(
        run.id,
        run.material_id,
        run.temperature[rows],
        run.temperature_unit,
        run.time === nothing ? nothing : run.time[rows],
        run.time_unit,
        [ThermalSignal(s.spec, s.values[rows]) for s in run.signals],
        [ThermalSegment(block.index, block.label, 1:length(rows))],
        copy(run.metadata),
        provenance,
    )
end

"""
    thermal_segment_diagnostics(run)

Summarize each block's row count, temperature range, temperature reversals, duration
and descriptive ordinary least-squares slope in K/min. Time-dependent fields are
`nothing` without recorded time; the slope is also `nothing` for a one-point block.
These checks do not classify program phases or establish kinetic applicability.
"""
function thermal_segment_diagnostics(run::ThermalRun)
    temperature = temperature_K(run)
    times = time_min(run)
    return map(run.segments) do block
        T = view(temperature, block.rows)
        ΔT = diff(T)
        t = times === nothing ? nothing : view(times, block.rows)
        duration = t === nothing ? nothing : last(t) - first(t)
        slope = if t === nothing || length(t) < 2
            nothing
        else
            centered_t = t .- mean(t)
            denominator = sum(abs2, centered_t)
            if denominator > 0 && isfinite(denominator)
                estimate = sum(centered_t .* (T .- mean(T))) / denominator
                isfinite(estimate) ? estimate : nothing
            else
                nothing
            end
        end
        return (;
            index=block.index,
            label=block.label,
            row_count=length(T),
            temperature_min_K=minimum(T),
            temperature_max_K=maximum(T),
            duration_min=duration,
            slope_K_per_min=slope,
            increasing_steps=count(>(0), ΔT),
            decreasing_steps=count(<(0), ΔT),
            flat_steps=count(iszero, ΔT),
        )
    end
end
