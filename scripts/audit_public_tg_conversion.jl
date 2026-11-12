using IsoconversionalAnalysis
using SHA
using Statistics
using TOML

module PublicTGImport
include("audit_public_thermal_import.jl")
end
using .PublicTGImport: public_input_path, require_check, thermal_mapping

const TG_ROOT = normpath(joinpath(@__DIR__, ".."))
const TG_IMPORT_MANIFEST = "config/public_thermal_import.toml"
const TG_PROTOCOL = "config/public_tg_conversion.toml"
const MIXCHAR_IMPORT_MANIFEST = "config/mixchar_thermal_import.toml"
const TG_REPORT = "docs/audits/public_tg_conversion.toml"

function window_rows(temperature_C, limits)
    lower, upper = limits
    indices = findall(T -> lower <= T <= upper, temperature_C)
    isempty(indices) && error("No recorded rows in [$lower, $upper] °C")
    rows = first(indices):last(indices)
    length(rows) == length(indices) || error(
        "[$lower, $upper] °C intersects separate blocks in acquisition order. Supply a row interval.",
    )
    return rows
end

function upward_crossings(curve, target)
    crossings = Int[]
    for i in 1:(length(curve.alpha) - 1)
        lower, upper = curve.alpha[i], curve.alpha[i + 1]
        if lower < target <= upper || (i == 1 && lower == target < upper)
            push!(crossings, i)
        end
    end
    result = Dict{String,Any}("alpha" => target, "count" => length(crossings))
    if !isempty(crossings)
        i = first(crossings)
        fraction = (target - curve.alpha[i]) / (curve.alpha[i + 1] - curve.alpha[i])
        left, right = curve.temperature_K[i], curve.temperature_K[i + 1]
        result["first_temperature_C"] = left + fraction * (right - left) - 273.15
        result["first_source_lines"] = curve.source_lines[i:(i + 1)]
    end
    return result
end

function reference_evidence(result)
    evidence = Dict{String,Any}("method" => String(result.method), "mass_mg" => result.mass)
    if result.segment_rows !== nothing
        evidence["source_lines"] = [first(result.source_lines), last(result.source_lines)]
        evidence["segment_rows"] = [first(result.segment_rows), last(result.segment_rows)]
        evidence["observed_mass_range_mg"] = collect(result.observed_mass_range)
    end
    return evidence
end

function case_evidence(curve, reference_curve, grid)
    crossings = [upward_crossings(curve, target) for target in grid]
    reference_crossings = [upward_crossings(reference_curve, target) for target in grid]
    shifts = Float64[]
    for (crossing, baseline) in zip(crossings, reference_crossings)
        if crossing["count"] == baseline["count"] == 1
            shift = crossing["first_temperature_C"] - baseline["first_temperature_C"]
            crossing["shift_from_plateau_K"] = shift
            push!(shifts, shift)
        end
    end
    a_min, a_max = extrema(curve.alpha)
    a_common = view(reference_curve.alpha, curve.reaction_rows)
    m0, mf = curve.initial_reference.mass, curve.terminal_reference.mass
    record = Dict{String,Any}(
        "rows" => length(curve.mass),
        "reaction_rows" => [first(curve.reaction_rows), last(curve.reaction_rows)],
        "reaction_source_lines" => [first(curve.source_lines), last(curve.source_lines)],
        "reaction_temperature_C" => collect(extrema(curve.temperature_K) .- 273.15),
        "initial_mass_mg" => m0,
        "terminal_mass_mg" => mf,
        "apparent_terminal_fraction" => mf / m0,
        "initial_reference" => reference_evidence(curve.initial_reference),
        "terminal_reference" => reference_evidence(curve.terminal_reference),
        "negative_terminal_reference" => mf < 0,
        "negative_mass_count" => curve.negative_mass_count,
        "alpha_first" => first(curve.alpha),
        "alpha_last" => last(curve.alpha),
        "alpha_min" => a_min,
        "alpha_max" => a_max,
        "decreasing_alpha_steps" => curve.decreasing_steps,
        "out_of_unit_interval_rows" => curve.outside_unit_interval_count,
        "max_abs_alpha_change_shared_rows" => maximum(abs, curve.alpha .- a_common),
        "unique_crossing_count" => count(c -> c["count"] == 1, crossings),
        "ambiguous_crossing_count" => count(c -> c["count"] > 1, crossings),
        "missing_crossing_count" => count(c -> c["count"] == 0, crossings),
        "crossings" => crossings,
    )
    isempty(shifts) || (record["max_abs_crossing_shift_K"] = maximum(abs, shifts))
    return record
end

function import_tg(expected, manifest)
    profile = manifest["profiles"][expected["profile"]]
    path = public_input_path(TG_ROOT, expected["path"])
    mapping = thermal_mapping(manifest, profile, Dict{String,String}())
    return read_thermal_delimited(
        path,
        mapping;
        id=expected["id"],
        material_id=expected["material"],
        delimiter=only(profile["delimiter"]),
        encoding=Symbol(profile["encoding"]),
        expected_sha256=expected["sha256"],
    )
end

function validate_case(curve, selected)
    id, rows = curve.run_id, curve.reaction_rows
    m0, mf = curve.initial_reference.mass, curve.terminal_reference.mass
    require_check(
        curve.source_lines == selected.provenance.source_lines[rows] &&
            curve.mass == signal_values(selected, curve.signal_id)[rows] &&
            curve.temperature_K == temperature_K(selected)[rows] &&
            curve.time_min == time_min(selected)[rows],
        id,
        "TG conversion altered observed rows or provenance",
    )
    reconstructed = m0 .- (m0 - mf) .* curve.alpha
    require_check(
        maximum(abs, reconstructed .- curve.mass) <= 1e-12, id, "mass closure failed"
    )
    return nothing
end

function tg_case_curves(run, expected, protocol)
    segment = expected["selected_segment"]
    selected = select_thermal_segment(run, segment)
    require_check(
        selected.temperature_unit == :C, run.id, "protocol windows require Celsius"
    )
    temperature = selected.temperature
    n = length(temperature)
    windows = protocol["windows"][expected["group"]]
    initial = window_rows(temperature, windows["initial_C"])
    terminal = window_rows(temperature, windows["terminal_C"])
    later_initial = window_rows(temperature, windows["later_initial_C"])
    earlier_terminal = window_rows(temperature, windows["earlier_terminal_C"])
    reaction = window_rows(temperature, windows["reaction_C"])
    convert(rows, m0, mf) = tg_conversion(
        run,
        :mass;
        segment,
        reaction_rows=rows,
        initial_reference=m0,
        terminal_reference=mf,
    )
    curves = Dict(
        "record_endpoints" => convert(1:n, 1:1, n:n),
        "plateau_medians" => convert(1:n, initial, terminal),
        "reaction_window" => convert(reaction, initial, terminal),
        "later_initial" => convert(1:n, later_initial, terminal),
        "earlier_terminal" => convert(1:n, initial, earlier_terminal),
        "window_and_references" => convert(reaction, later_initial, earlier_terminal),
    )
    if expected["group"] == "experimental_tg"
        curves["zero_terminal"] = convert(1:n, initial, 0.0)
    end
    reference = curves[protocol["comparison_case"]]
    require_check(reference.mass_unit == :mg, run.id, "reference units must be mg")
    require_check(
        curves["reaction_window"].alpha == reference.alpha[reaction],
        run.id,
        "changing the interval without changing references altered overlapping conversion",
    )
    require_check(
        first(curves["record_endpoints"].alpha) == 0 &&
            last(curves["record_endpoints"].alpha) == 1,
        run.id,
        "endpoint normalization failed to recover endpoints",
    )
    for curve in values(curves)
        validate_case(curve, selected)
    end
    return (; selected, curves, reference)
end

function audit_tg_run(expected, manifest, protocol)
    run = import_tg(expected, manifest)
    selected, curves, reference = tg_case_curves(run, expected, protocol)
    measured_rate = only(thermal_segment_diagnostics(selected)).slope_K_per_min
    return Dict{String,Any}(
        "id" => run.id,
        "material_id" => run.material_id,
        "group" => expected["group"],
        "source_path" => expected["path"],
        "source_sha256" => run.provenance.sha256,
        "segment_index" => only(selected.segments).index,
        "rows" => length(selected.temperature),
        "measured_rate_K_per_min" => measured_rate,
        "cases" => Dict(
            name => case_evidence(curve, reference, protocol["alpha_grid"]) for
            (name, curve) in curves
        ),
    )
end

function audit_mixchar_losses(config, protocol)
    temperature = ThermalColumnSpec(
        config["temperature_column"], Symbol(config["temperature_unit"])
    )
    signal = ThermalSignalSpec(
        :mass_loss,
        ThermalColumnSpec(config["signal_column"], Symbol(config["signal_unit"])),
        Symbol(config["signal_quantity"]),
    )
    mapping = ThermalColumnMapping(temperature, [signal])
    return map(config["runs"]) do expected
        id = expected["id"]
        rda_path = public_input_path(TG_ROOT, expected["rda_path"])
        rda_hash = bytes2hex(open(sha256, rda_path))
        require_check(rda_hash == expected["rda_sha256"], id, "unexpected RDA source hash")
        csv_path = joinpath(TG_ROOT, expected["csv_path"])
        isfile(csv_path) ||
            error("Missing $csv_path. Run scripts/prepare_mixchar_examples.py.")
        run = read_thermal_delimited(
            csv_path,
            mapping;
            id,
            material_id=expected["material_id"],
            expected_sha256=expected["csv_sha256"],
        )
        require_check(
            length(run.temperature) == expected["rows"], id, "unexpected row count"
        )
        require_check(time_min(run) === nothing, id, "unexpected recorded time")
        original_mass = expected["vignette_initial_mass_mg"]
        loss = signal_values(run, :mass_loss)
        fraction = -loss ./ original_mass
        require_check(
            all(x -> 0 <= x <= 1, fraction), id, "invalid original-sample loss fraction"
        )
        late = window_rows(run.temperature, protocol["windows"]["mixchar"]["late_C"])
        return Dict{String,Any}(
            "id" => id,
            "material_id" => run.material_id,
            "rda_path" => expected["rda_path"],
            "rda_sha256" => rda_hash,
            "csv_path" => expected["csv_path"],
            "csv_sha256" => run.provenance.sha256,
            "rows" => length(fraction),
            "recorded_time" => false,
            "reaction_conversion_defined" => false,
            "vignette_initial_mass_mg" => original_mass,
            "terminal_observed_temperature_C" => last(run.temperature),
            "terminal_observed_signed_mass_loss_mg" => last(loss),
            "terminal_original_sample_loss_fraction" => last(fraction),
            "late_window_median_original_sample_loss_fraction" =>
                median(view(fraction, late)),
            "late_window_source_lines" =>
                [run.provenance.source_lines[i] for i in (first(late), last(late))],
        )
    end
end

function common_crossings(runs, grid, name)
    return [
        a for (i, a) in enumerate(grid) if
        all(run -> run["cases"][name]["crossings"][i]["count"] == 1, runs)
    ]
end

function material_summary(runs, grid)
    return map(sort(unique(run["material_id"] for run in runs))) do material
        group = filter(run -> run["material_id"] == material, runs)
        return Dict{String,Any}(
            "material_id" => material,
            "runs" => length(group),
            "unique_crossings_in_all_runs" => Dict(
                name => common_crossings(group, grid, name) for
                name in keys(first(group)["cases"])
            ),
        )
    end
end

function main()
    picnik_only = "--picnik-only" in ARGS
    paths = filter(!=("--picnik-only"), ARGS)
    length(paths) <= 1 || error(
        "Usage: julia --project=. scripts/audit_public_tg_conversion.jl [--picnik-only] [report.toml]",
    )
    default_report = picnik_only ? "results/public_tg_conversion.toml" : TG_REPORT
    report_path = isempty(paths) ? joinpath(TG_ROOT, default_report) : abspath(only(paths))
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    protocol = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    protocol["schema_version"] == 1 || error("Unsupported TG protocol schema")
    groups = ("experimental_tg", "synthetic_tg")
    runs = [
        audit_tg_run(r, manifest, protocol) for
        r in manifest["runs"] if r["group"] in groups
    ]
    count(r -> r["group"] == "synthetic_tg", runs) == 13 ||
        error("Expected thirteen synthetic TG runs")
    count(r -> r["group"] == "experimental_tg", runs) == 5 ||
        error("Expected five experimental TG runs")
    mixchar = if picnik_only
        Dict{String,Any}[]
    else
        config = TOML.parsefile(joinpath(TG_ROOT, MIXCHAR_IMPORT_MANIFEST))
        losses = audit_mixchar_losses(config, protocol)
        length(losses) == 2 || error("Expected two mixchar runs")
        losses
    end
    evidence = Dict{String,Any}(
        "schema_version" => 1,
        "import_manifest" => TG_IMPORT_MANIFEST,
        "protocol" => TG_PROTOCOL,
        "alpha_grid" => protocol["alpha_grid"],
        "comparison_case" => protocol["comparison_case"],
        "tg_conversion_runs" => length(runs),
        "mixchar_original_sample_loss_runs" => length(mixchar),
        "tg_conversion_cases" => sum(r -> length(r["cases"]), runs),
        "material_summary" => material_summary(runs, protocol["alpha_grid"]),
        "runs" => runs,
        "mixchar" => mixchar,
    )
    picnik_only || (evidence["mixchar_manifest"] = MIXCHAR_IMPORT_MANIFEST)
    mkpath(dirname(report_path))
    open(report_path, "w") do io
        return TOML.print(io, evidence; sorted=true)
    end
    println(
        "Verified $(evidence["tg_conversion_cases"]) TG cases across $(length(runs)) runs and $(length(mixchar)) mixchar loss traces.",
    )
    return println("Wrote $report_path")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
