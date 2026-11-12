using IsoconversionalAnalysis
using Statistics
using TOML

module PublicTGConversionAudit
include("audit_public_tg_conversion.jl")
end

using .PublicTGConversionAudit:
    TG_ROOT, TG_IMPORT_MANIFEST, TG_PROTOCOL, import_tg, tg_case_curves, require_check

const SUPPORT_REPORT = "docs/audits/public_tg_conversion_support.toml"
const GROUP_RUN_COUNTS = Dict(
    "Constant_E" => 5, "Two_Steps" => 4, "Variable_E" => 4, "polypropylene" => 5
)
const PP_NEGATIVE_TERMINAL_RUN = "picnik_PP_10_B3"
const CENTRAL_ALPHA_BOUNDS = (0.10, 0.90)

function program_diagnostics(curve, rows)
    t = view(curve.time_min, rows)
    T = view(curve.temperature_K, rows)
    dt = diff(t)
    dT = diff(T)
    center_t = t .- mean(t)
    center_T = T .- mean(T)
    slope = sum(center_t .* center_T) / sum(abs2, center_t)
    residual = center_T .- slope .* center_t
    return Dict{String,Any}(
        "rows" => length(rows),
        "source_line_range" =>
            [curve.source_lines[first(rows)], curve.source_lines[last(rows)]],
        "temperature_range_K" => collect(extrema(T)),
        "elapsed_min" => last(t) - first(t),
        "OLS_slope_K_per_min" => slope,
        "OLS_temperature_residual_RMSE_K" => sqrt(mean(abs2, residual)),
        "OLS_max_abs_temperature_residual_K" => maximum(abs, residual),
        "step_slope_range_K_per_min" => collect(extrema(dT ./ dt)),
        "falling_temperature_steps" => count(<(0), dT),
        "flat_temperature_steps" => count(==(0), dT),
    )
end

function central_program_diagnostics(curve)
    diagnostics = [tg_conversion_crossings(curve, alpha) for alpha in CENTRAL_ALPHA_BOUNDS]
    if !all(d -> d.status === :unique_upward, diagnostics)
        return Dict("boundary_statuses" => [String(d.status) for d in diagnostics])
    end
    lower, upper = (only(d.crossings).segment_rows for d in diagnostics)
    offset = first(curve.reaction_rows) - 1
    rows = (first(lower) - offset):(last(upper) - offset)
    result = program_diagnostics(curve, rows)
    result["boundary_statuses"] = [String(d.status) for d in diagnostics]
    midpoint = tg_conversion_crossings(curve, 0.50)
    if midpoint.status === :unique_upward
        crossings = [
            only(d.crossings) for d in (first(diagnostics), midpoint, last(diagnostics))
        ]
        t = [event.time_min::Float64 for event in crossings]
        T = [event.temperature_K::Float64 for event in crossings]
        result["crossing_secant_alpha_nodes"] = [
            first(CENTRAL_ALPHA_BOUNDS), 0.50, last(CENTRAL_ALPHA_BOUNDS)
        ]
        result["crossing_secant_K_per_min"] = diff(T) ./ diff(t)
    end
    return result
end

function load_run(spec, manifest, protocol)
    run = import_tg(spec, manifest)
    selected, curves, reference = tg_case_curves(run, spec, protocol)
    require_check(
        reference.time_min !== nothing &&
            reference.reaction_rows == 1:length(selected.temperature),
        run.id,
        "program diagnostics require the complete recorded time history",
    )
    full_program = program_diagnostics(reference, eachindex(reference.temperature_K))
    slope = only(thermal_segment_diagnostics(selected)).slope_K_per_min
    require_check(
        isapprox(full_program["OLS_slope_K_per_min"], slope; rtol=1e-12),
        run.id,
        "full-history slope disagrees with the import audit",
    )
    evidence = Dict{String,Any}(
        "id" => run.id,
        "source_path" => spec["path"],
        "source_sha256" => run.provenance.sha256,
        "segment_index" => only(selected.segments).index,
        "comparison_terminal_reference_mg" => reference.terminal_reference.mass,
        "full_program" => full_program,
        "central_alpha_program" => central_program_diagnostics(reference),
    )
    return (;
        id=run.id, material=run.material_id, group=spec["group"], slope, curves, evidence
    )
end

function event_evidence(event, curve)
    rows = event.segment_rows
    offset = first(curve.reaction_rows) - 1
    left, right = first(rows) - offset, last(rows) - offset
    require_check(
        event.source_lines == (curve.source_lines[left], curve.source_lines[right]),
        curve.run_id,
        "crossing provenance disagrees with the observed source lines",
    )
    evidence = Dict{String,Any}(
        "kind" => String(event.kind),
        "segment_row_range" => [first(rows), last(rows)],
        "source_line_range" => collect(event.source_lines),
        "temperature_range_K" => collect(event.temperature_range_K),
    )
    if event.temperature_K !== nothing
        evidence["temperature_K"] = event.temperature_K
    end
    if event.time_min !== nothing
        evidence["time_min"] = event.time_min
        if right > left
            evidence["observed_bracket_slope_K_per_min"] =
                (curve.temperature_K[right] - curve.temperature_K[left]) /
                (curve.time_min[right] - curve.time_min[left])
        end
    end
    return evidence
end

function coverage_point(point, run_ids)
    upward = [
        id for (id, diagnostic) in zip(run_ids, point.crossings) if
        diagnostic.status === :unique_upward
    ]
    require_check(
        length(upward) == point.upward_run_count,
        point.alpha,
        "per-run upward support count changed",
    )
    evidence = Dict{String,Any}(
        "alpha" => point.alpha,
        "common_upward" => point.common_upward,
        "upward_run_count" => length(upward),
        "upward_run_ids" => upward,
    )
    if point.common_upward
        T = [only(d.crossings).temperature_K::Float64 for d in point.crossings]
        cold, hot = argmin(T), argmax(T)
        require_check(
            point.temperature_min_K == T[cold] && point.temperature_max_K == T[hot],
            point.alpha,
            "temperature extrema disagree with the crossing events",
        )
        evidence["temperature_min_K"] = point.temperature_min_K
        evidence["temperature_mean_K"] = point.temperature_mean_K
        evidence["temperature_max_K"] = point.temperature_max_K
        evidence["temperature_span_K"] = point.temperature_span_K
        evidence["temperature_min_run_id"] = run_ids[cold]
        evidence["temperature_max_run_id"] = run_ids[hot]
        evidence["reciprocal_temperature_span_K_inv"] =
            point.reciprocal_temperature_span_K_inv
        evidence["reciprocal_temperature_Sxx_K_inv_sq"] =
            point.reciprocal_temperature_sxx_K_inv2
    end
    return evidence
end

function add_reference_comparison!(points, reference; paired=false)
    for (point, base) in zip(points, reference)
        if point["common_upward"] && base["common_upward"]
            point["temperature_span_change_K"] =
                point["temperature_span_K"] - base["temperature_span_K"]
            if paired
                temperatures(p) =
                    [only(row["events"])["temperature_K"] for row in p["crossings"]]
                point["max_abs_paired_crossing_shift_K"] = maximum(
                    abs, temperatures(point) .- temperatures(base)
                )
            end
        end
    end
    return nothing
end

common_alpha(points) = [p["alpha"] for p in points if p["common_upward"]]

function case_evidence(runs, case, targets, reference_case)
    curves = [run.curves[case] for run in runs]
    support = tg_conversion_support(curves, targets)
    require_check(support.run_ids == [run.id for run in runs], case, "run order changed")
    require_check(
        isempty(support.missing_time_run_ids), case, "public runs require recorded time"
    )
    points = map(support.points) do point
        evidence = coverage_point(point, support.run_ids)
        evidence["crossings"] = [
            Dict{String,Any}(
                "run_id" => curve.run_id,
                "status" => String(diagnostic.status),
                "events" =>
                    [event_evidence(event, curve) for event in diagnostic.crossings],
            ) for (curve, diagnostic) in zip(curves, point.crossings)
        ]
        return evidence
    end
    omitted = if case == reference_case
        support.run_ids
    elseif PP_NEGATIVE_TERMINAL_RUN in support.run_ids
        [PP_NEGATIVE_TERMINAL_RUN]
    else
        String[]
    end
    subsets = map(omitted) do id
        selected = filter(c -> c.run_id != id, curves)
        subset = tg_conversion_support(selected, targets)
        subset_points = [coverage_point(point, subset.run_ids) for point in subset.points]
        add_reference_comparison!(subset_points, points)
        return Dict{String,Any}(
            "id" => "omit_" * id,
            "omitted_run_id" => id,
            "run_ids" => subset.run_ids,
            "common_alpha" => common_alpha(subset_points),
            "added_vs_full_alpha" =>
                setdiff(common_alpha(subset_points), common_alpha(points)),
            "points" => subset_points,
        )
    end
    return Dict{String,Any}(
        "id" => case,
        "common_alpha" => common_alpha(points),
        "points" => points,
        "subsets" => subsets,
    )
end

function group_evidence(material, runs, targets, reference_case)
    selected = sort([run for run in runs if run.material == material]; by=run -> run.slope)
    require_check(
        length(selected) == GROUP_RUN_COUNTS[material], material, "run roster changed"
    )
    ids = [run.id for run in selected]
    case_ids = sort!(collect(keys(first(selected).curves)))
    require_check(
        all(run -> Set(keys(run.curves)) == Set(case_ids), selected),
        material,
        "conversion scenarios differ across runs",
    )
    cases = [case_evidence(selected, id, targets, reference_case) for id in case_ids]
    reference = only(filter(case -> case["id"] == reference_case, cases))
    for case in cases
        add_reference_comparison!(case["points"], reference["points"]; paired=true)
        case["shared_with_reference_alpha"] = intersect(
            case["common_alpha"], reference["common_alpha"]
        )
        case["added_vs_reference_alpha"] = setdiff(
            case["common_alpha"], reference["common_alpha"]
        )
        case["lost_vs_reference_alpha"] = setdiff(
            reference["common_alpha"], case["common_alpha"]
        )
    end
    common_omissions = reduce(
        intersect, [[s["omitted_run_id"] for s in case["subsets"]] for case in cases]
    )
    omission_intersections = map(common_omissions) do id
        subsets = [
            only(filter(s -> s["omitted_run_id"] == id, case["subsets"])) for case in cases
        ]
        return Dict{String,Any}(
            "omitted_run_id" => id,
            "common_across_all_cases_alpha" =>
                reduce(intersect, [s["common_alpha"] for s in subsets]),
        )
    end
    return Dict{String,Any}(
        "material_id" => material,
        "dataset_group" => only(unique(run.group for run in selected)),
        "run_ids" => ids,
        "runs" => [run.evidence for run in selected],
        "cases" => cases,
        "common_across_all_cases_alpha" =>
            reduce(intersect, [case["common_alpha"] for case in cases]),
        "omission_intersections" => omission_intersections,
    )
end

function main()
    length(ARGS) <= 1 || error(
        "Usage: julia --project=. scripts/audit_public_tg_conversion_support.jl [report.toml]",
    )
    report = isempty(ARGS) ? joinpath(TG_ROOT, SUPPORT_REPORT) : abspath(only(ARGS))
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    protocol = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    require_check(
        protocol["schema_version"] == 1, "protocol", "unexpected TG protocol schema"
    )
    targets = protocol["alpha_grid"]
    reference_case = protocol["comparison_case"]
    runs = [
        load_run(spec, manifest, protocol) for
        spec in manifest["runs"] if spec["group"] in ("synthetic_tg", "experimental_tg")
    ]
    materials = sort!(collect(keys(GROUP_RUN_COUNTS)))
    require_check(
        Set(run.material for run in runs) == Set(materials), "roster", "materials changed"
    )
    pp10 = only(filter(run -> run.id == PP_NEGATIVE_TERMINAL_RUN, runs))
    require_check(
        pp10.curves[reference_case].terminal_reference.mass < 0,
        pp10.id,
        "expected a negative apparent terminal reference in the sensitivity trace",
    )
    groups = [
        group_evidence(material, runs, targets, reference_case) for material in materials
    ]
    evidence = Dict{String,Any}(
        "schema_version" => 1,
        "raw_import_manifest" => TG_IMPORT_MANIFEST,
        "tg_protocol" => TG_PROTOCOL,
        "alpha_grid" => targets,
        "comparison_case" => reference_case,
        "central_program_alpha_bounds" => collect(CENTRAL_ALPHA_BOUNDS),
        "run_count" => length(runs),
        "group_count" => length(groups),
        "group_case_count" => sum(length(group["cases"]) for group in groups),
        "subset_case_count" =>
            sum(length(case["subsets"]) for group in groups for case in group["cases"]),
        "full_group_crossing_evaluations" =>
            length(targets) * sum(length(run.curves) for run in runs),
        "groups" => groups,
    )
    open(report, "w") do io
        return TOML.print(io, evidence; sorted=true)
    end
    println(
        "Verified crossing support for $(length(runs)) TG runs in $(length(groups)) groups."
    )
    return println("Wrote $report")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
