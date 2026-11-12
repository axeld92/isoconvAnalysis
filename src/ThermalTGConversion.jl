"""
A mass reference used by [`tg_conversion`](@ref) in the selected signal's native unit.
`method` records `:fixed` for a supplied mass or `:median` for an observed window.
Observed windows include segment-local row indices, source lines, and observed mass range.
"""
struct TGReferenceEstimate
    mass::Float64
    method::Symbol
    segment_rows::Union{Nothing,UnitRange{Int}}
    source_lines::Vector{Int}
    observed_mass_range::Union{Nothing,Tuple{Float64,Float64}}
end

"""
TG reaction progress conditional on a mass-loss interval and two mass references.
`reaction_rows` indexes the chosen acquisition block. `mass` is the measured signal in
`mass_unit`. `temperature_K` and optional `time_min` preserve measurement order and the
clock origin. `alpha` is dimensionless and retains excursions outside zero to one.
`decreasing_steps` and `outside_unit_interval_count` diagnose the resulting curve.
`negative_mass_count` flags negative readings in the selected interval.
"""
struct TGConversionResult
    run_id::String
    material_id::Union{Nothing,String}
    segment_index::Int
    segment_label::Union{Nothing,String}
    signal_id::Symbol
    mass_unit::Symbol
    reaction_rows::UnitRange{Int}
    source_path::String
    source_sha256::String
    source_lines::Vector{Int}
    temperature_K::Vector{Float64}
    time_min::Union{Nothing,Vector{Float64}}
    mass::Vector{Float64}
    initial_reference::TGReferenceEstimate
    terminal_reference::TGReferenceEstimate
    alpha::Vector{Float64}
    decreasing_steps::Int
    outside_unit_interval_count::Int
    negative_mass_count::Int
end

function _check_tg_rows(rows, segment_length, label)
    !isempty(rows) && 1 <= first(rows) <= last(rows) <= segment_length ||
        throw(ArgumentError("$label rows must be inside the selected segment"))
    return nothing
end

function _tg_reference(values, lines, reference::Real, label)
    mass = Float64(reference)
    isfinite(mass) || throw(ArgumentError("$label reference mass must be finite"))
    return TGReferenceEstimate(mass, :fixed, nothing, Int[], nothing)
end

function _tg_reference(values, lines, rows::UnitRange{Int}, label)
    _check_tg_rows(rows, length(values), label)
    observed = view(values, rows)
    all(isfinite, observed) ||
        throw(ArgumentError("$label reference contains nonfinite mass"))
    return TGReferenceEstimate(
        median(observed), :median, rows, lines[rows], extrema(observed)
    )
end

"""
    tg_conversion(run, signal_id; segment, reaction_rows, initial_reference, terminal_reference)

Compute `alpha = (initial_mass - mass) / (initial_mass - terminal_mass)` for a remaining
mass channel (`:mass`). Choose the original segment ordinal and a contiguous inclusive
`reaction_rows` range indexed within that segment. Supply each reference as a finite mass
in the channel's native unit or as an inclusive segment-local row range whose median mass
defines the reference. A one-row range uses the recorded endpoint. When both references
use observed windows, the initial window must precede the terminal window.

The initial reference must be positive and exceed the terminal reference. A negative
measured terminal mass can define apparent mass loss and warrants checking the balance
baseline before interpreting physical residue. The calculation preserves acquisition order,
unclipped progress, measured temperature reversals, and optional recorded time. Interpreting
the result as chemical conversion or reaction completion requires separate evidence.
"""
function tg_conversion(
    run::ThermalRun,
    signal_id::Symbol;
    segment::Integer,
    reaction_rows::UnitRange{Int},
    initial_reference::Union{Real,UnitRange{Int}},
    terminal_reference::Union{Real,UnitRange{Int}},
)
    selected = select_thermal_segment(run, segment)
    signal_index = findfirst(signal -> signal.spec.id === signal_id, selected.signals)
    signal_index === nothing &&
        throw(ArgumentError("run $(run.id) has no signal $signal_id"))
    signal = selected.signals[signal_index]
    signal.spec.quantity === :mass ||
        throw(ArgumentError("TG conversion requires a remaining mass channel tagged :mass"))
    _check_tg_rows(reaction_rows, length(selected.temperature), "reaction")
    if initial_reference isa UnitRange && terminal_reference isa UnitRange
        last(initial_reference) < first(terminal_reference) || throw(
            ArgumentError("initial reference rows must precede terminal reference rows")
        )
    end
    lines = selected.provenance.source_lines
    initial = _tg_reference(signal.values, lines, initial_reference, "initial")
    terminal = _tg_reference(signal.values, lines, terminal_reference, "terminal")
    isfinite(initial.mass) &&
    isfinite(terminal.mass) &&
    initial.mass > max(terminal.mass, 0) || throw(
        ArgumentError(
            "TG references require positive initial mass above finite terminal mass"
        ),
    )
    denominator = initial.mass - terminal.mass
    isfinite(denominator) || throw(ArgumentError("TG mass difference must be finite"))
    mass = signal.values[reaction_rows]
    alpha = (initial.mass .- mass) ./ denominator
    all(isfinite, alpha) || throw(ArgumentError("TG conversion contains nonfinite values"))
    t = time_min(selected)
    block = only(selected.segments)
    return TGConversionResult(
        selected.id,
        selected.material_id,
        block.index,
        block.label,
        signal_id,
        signal.spec.column.unit,
        reaction_rows,
        selected.provenance.source_path,
        selected.provenance.sha256,
        lines[reaction_rows],
        temperature_K(selected)[reaction_rows],
        t === nothing ? nothing : t[reaction_rows],
        mass,
        initial,
        terminal,
        alpha,
        count(<(0), diff(alpha)),
        count(a -> !(0 <= a <= 1), alpha),
        count(<(0), mass),
    )
end
