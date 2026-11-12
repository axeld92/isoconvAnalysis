using IsoconversionalAnalysis
using TOML

module PublicVyazovkinConversionAudit
include("audit_public_tg_conversion.jl")
end

using .PublicVyazovkinConversionAudit:
    TG_ROOT, TG_IMPORT_MANIFEST, TG_PROTOCOL, import_tg, tg_case_curves, require_check

const VYAZOVKIN_AUDIT_PROTOCOL = "config/public_tg_vyazovkin.toml"
const VYAZOVKIN_AUDIT_REPORT = "docs/audits/public_tg_vyazovkin.toml"
const TG_MATERIAL_COUNTS = Dict(
    "Constant_E" => 5, "Two_Steps" => 4, "Variable_E" => 4, "polypropylene" => 5
)

function recorded_path_evidence!(evidence, curve, interval)
    t = interval.time_min
    isempty(t) && return evidence
    lines = interval.source_lines
    raw_nodes = 0
    for (i, (left, right)) in enumerate(lines)
        left == right || continue
        row = searchsortedfirst(curve.source_lines, left)
        require_check(
            curve.source_lines[row] == left &&
                curve.time_min[row] == t[i] &&
                curve.temperature_K[row] == interval.temperature_K[i] &&
                curve.alpha[row] == interval.alpha[i],
            curve.run_id,
            "a measured node changed inside the clipped interval",
        )
        raw_nodes += 1
    end
    evidence["boundary_source_lines"] = [collect(first(lines)), collect(last(lines))]
    evidence["boundary_times_min"] = [first(t), last(t)]
    evidence["boundary_temperatures_K"] = [
        first(interval.temperature_K), last(interval.temperature_K)
    ]
    evidence["duration_min"] = last(t) - first(t)
    evidence["clipped_node_count"] = length(t)
    evidence["measured_node_count"] = raw_nodes
    first_row = searchsortedfirst(curve.source_lines, first(lines)[1])
    last_row = searchsortedfirst(curve.source_lines, last(lines)[2])
    evidence["single_acquisition_bracket"] = last_row - first_row == 1
    evidence["max_time_step_min"] = maximum(diff(t))
    evidence["temperature_range_K"] = collect(extrema(interval.temperature_K))
    evidence["falling_temperature_steps"] = count(<(0), diff(interval.temperature_K))
    evidence["reversal_steps"] = interval.reversal_steps
    evidence["plateau_steps"] = interval.plateau_steps
    evidence["max_backtracking_alpha"] = interval.max_backtracking_alpha
    evidence["plateau_duration_min"] = interval.plateau_duration_min
    return evidence
end

function quadrature_evidence!(evidence, interval, options)
    probes = [
        tg_recorded_time_quadrature(
            interval,
            probe;
            legacy_alpha_nodes=options["legacy_alpha_nodes"],
            subdivisions_per_step=options["subdivisions_per_step"],
        ) for probe in options["probe_energies_kJ_per_mol"]
    ]
    native = getfield.(probes, :recorded_log_integral)
    legacy = getfield.(probes, :legacy_log_integral)
    refined = getfield.(probes, :refined_log_integral)
    doubled = getfield.(probes, :doubled_log_integral)
    evidence["log_recorded"] = native
    evidence["log_legacy_uniform_alpha"] = legacy
    evidence["log_refined_time"] = refined
    evidence["log_doubled_time"] = doubled
    evidence["recorded_relative_difference"] = expm1.(native .- doubled)
    evidence["legacy_relative_difference"] = expm1.(legacy .- doubled)
    evidence["refinement_relative_difference"] = expm1.(refined .- doubled)
    return evidence
end

function interval_evidence(curve, alpha, delta, options)
    interval = tg_recorded_time_interval(curve, alpha; delta_alpha=delta)
    evidence = Dict{String,Any}(
        "alpha" => alpha,
        "alpha_lower" => interval.alpha_lower,
        "status" => String(interval.status),
        "lower_crossing_status" => String(interval.lower_crossing_status),
        "upper_crossing_status" => String(interval.upper_crossing_status),
    )
    recorded_path_evidence!(evidence, curve, interval)
    if interval.status === :supported
        quadrature_evidence!(evidence, interval, options)
    end
    return evidence
end

function load_vyazovkin_run(spec, manifest, protocol, options)
    run = import_tg(spec, manifest)
    _, curves, reference = tg_case_curves(run, spec, protocol)
    cases = map(sort!(collect(keys(curves)))) do case_id
        curve = curves[case_id]
        windows = map(options["delta_alpha_grid"]) do delta
            return Dict{String,Any}(
                "delta_alpha" => delta,
                "points" => [
                    interval_evidence(curve, alpha, delta, options) for
                    alpha in protocol["alpha_grid"]
                ],
            )
        end
        return Dict("id" => case_id, "windows" => windows)
    end
    return Dict{String,Any}(
        "run_id" => run.id,
        "material_id" => run.material_id,
        "source_sha256" => run.provenance.sha256,
        "source_path" => spec["path"],
        "segment_index" => reference.segment_index,
        "reference_terminal_mass_mg" => reference.terminal_reference.mass,
        "cases" => cases,
    )
end

supported(point) = point["status"] == "supported"

function common_targets(points, targets)
    return [targets[i] for i in eachindex(targets) if all(run -> supported(run[i]), points)]
end

function status_counts(points)
    counts = Dict{String,Int}()
    for point in Iterators.flatten(points)
        status = point["status"]
        counts[status] = get(counts, status, 0) + 1
    end
    return counts
end

function group_window_summary(points, delta, targets, options, run_ids)
    common = common_targets(points, targets)
    lower, upper = options["central_alpha_bounds"]
    active = filter(supported, collect(Iterators.flatten(points)))
    evidence = Dict{String,Any}(
        "delta_alpha" => delta,
        "status_counts" => status_counts(points),
        "common_supported_alpha" => common,
        "central_supported_alpha" => filter(a -> lower <= a <= upper, common),
    )
    if !isempty(active)
        evidence["min_measured_node_count"] = minimum(
            point -> point["measured_node_count"], active
        )
        evidence["single_bracket_supported_count"] = count(
            point -> point["single_acquisition_bracket"], active
        )
        evidence["supported_duration_range_min"] = collect(
            extrema(point -> point["duration_min"], active)
        )
        for (label, key) in (
            ("recorded", "recorded_relative_difference"),
            ("legacy", "legacy_relative_difference"),
            ("refinement", "refinement_relative_difference"),
        )
            evidence["max_abs_$(label)_relative_difference"] = [
                maximum(point -> abs(point[key][i]), active) for
                i in eachindex(options["probe_energies_kJ_per_mol"])
            ]
        end
    end
    omitted = options["polypropylene_omitted_run"]
    if omitted in run_ids
        retained = [i for i in eachindex(run_ids) if run_ids[i] != omitted]
        evidence["omit_pp10_common_supported_alpha"] = common_targets(
            points[retained], targets
        )
    end
    return evidence
end

function group_summary(material, runs, targets, options)
    selected = filter(run -> run["material_id"] == material, runs)
    require_check(
        length(selected) == TG_MATERIAL_COUNTS[material], material, "unexpected run roster"
    )
    run_ids = [run["run_id"] for run in selected]
    cases = map(first(selected)["cases"]) do first_case
        case_id = first_case["id"]
        per_run = [
            only(filter(case -> case["id"] == case_id, run["cases"]))["windows"] for
            run in selected
        ]
        windows = map(enumerate(options["delta_alpha_grid"])) do (i, delta)
            points = [window[i]["points"] for window in per_run]
            return group_window_summary(points, delta, targets, options, run_ids)
        end
        return Dict("id" => case_id, "windows" => windows)
    end
    return Dict("material_id" => material, "run_ids" => run_ids, "cases" => cases)
end

function main()
    length(ARGS) <= 1 ||
        error("Usage: julia --project=. scripts/audit_public_tg_vyazovkin.jl [report.toml]")
    report = isempty(ARGS) ? joinpath(TG_ROOT, VYAZOVKIN_AUDIT_REPORT) : abspath(only(ARGS))
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    protocol = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    options = TOML.parsefile(joinpath(TG_ROOT, VYAZOVKIN_AUDIT_PROTOCOL))
    require_check(options["schema_version"] == 1, "protocol", "unknown audit protocol")
    specs = [
        spec for
        spec in manifest["runs"] if spec["group"] in ("synthetic_tg", "experimental_tg")
    ]
    runs = [load_vyazovkin_run(spec, manifest, protocol, options) for spec in specs]
    require_check(
        Set(run["material_id"] for run in runs) == Set(keys(TG_MATERIAL_COUNTS)),
        "roster",
        "unexpected materials",
    )
    omitted = only(
        filter(run -> run["run_id"] == options["polypropylene_omitted_run"], runs)
    )
    require_check(
        omitted["reference_terminal_mass_mg"] < 0,
        omitted["run_id"],
        "expected a negative apparent terminal reference in the omission trace",
    )
    targets = protocol["alpha_grid"]
    groups = [
        group_summary(material, runs, targets, options) for
        material in sort!(collect(keys(TG_MATERIAL_COUNTS)))
    ]
    windows = [
        window for run in runs for case in run["cases"] for window in case["windows"]
    ]
    evidence = Dict{String,Any}(
        "schema_version" => 1,
        "raw_import_manifest" => TG_IMPORT_MANIFEST,
        "tg_protocol" => TG_PROTOCOL,
        "vyazovkin_protocol" => VYAZOVKIN_AUDIT_PROTOCOL,
        "alpha_grid" => targets,
        "delta_alpha_grid" => options["delta_alpha_grid"],
        "probe_energies_kJ_per_mol" => options["probe_energies_kJ_per_mol"],
        "legacy_alpha_nodes" => options["legacy_alpha_nodes"],
        "subdivisions_per_step" => options["subdivisions_per_step"],
        "central_alpha_bounds" => options["central_alpha_bounds"],
        "polypropylene_omitted_run" => options["polypropylene_omitted_run"],
        "comparison_case" => protocol["comparison_case"],
        "run_count" => length(runs),
        "run_case_count" => sum(run -> length(run["cases"]), runs),
        "interval_count" => sum(window -> length(window["points"]), windows),
        "status_counts" => status_counts([window["points"] for window in windows]),
        "groups" => groups,
        "runs" => runs,
    )
    open(report, "w") do io
        return TOML.print(io, evidence; sorted=true)
    end
    println(
        "Checked $(evidence["interval_count"]) measured-time intervals on $(length(runs)) TG runs.",
    )
    return println("Wrote $report")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
