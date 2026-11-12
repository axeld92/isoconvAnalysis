using IsoconversionalAnalysis
using TOML

module PublicTGFriedmanInputs
include("audit_public_tg_conversion.jl")
end

using .PublicTGFriedmanInputs:
    TG_ROOT, TG_IMPORT_MANIFEST, TG_PROTOCOL, import_tg, tg_case_curves

const FRIEDMAN_PROTOCOL = "config/public_tg_friedman.toml"
const FRIEDMAN_REPORT = "docs/audits/public_tg_friedman.toml"

function friedman_clock_scale(curve, bounds)
    crossings = tg_conversion_crossings.(Ref(curve), bounds)
    all(c -> c.status === :unique_upward, crossings) ||
        error("$(curve.run_id) lacks an upward crossing at the derivative clock anchors")
    t = [only(c.crossings).time_min for c in crossings]
    all(!isnothing, t) && t[2] > t[1] || error(
        "$(curve.run_id) needs ordered recorded times at the derivative clock anchors"
    )
    return t[2] - t[1]
end

function friedman_case(curve, grid, settings, windows, minimum_points)
    return map(grid) do alpha
        return [
            tg_friedman_rate(
                curve,
                alpha;
                half_window_min=width,
                degree=setting["degree"],
                minimum_points_per_side=minimum_points,
            ) for (setting, width) in zip(settings, windows)
        ]
    end
end

function friedman_run(spec, manifest, conversion_protocol, protocol)
    run = import_tg(spec, manifest)
    _, cases, reference = tg_case_curves(run, spec, conversion_protocol)
    clock_scale = friedman_clock_scale(reference, protocol["clock_scale_alpha"])
    settings = protocol["settings"]
    windows = [clock_scale * setting["half_window_fraction"] for setting in settings]
    cases = Dict(
        id => friedman_case(
            curve,
            conversion_protocol["alpha_grid"],
            settings,
            windows,
            protocol["minimum_points_per_side"],
        ) for (id, curve) in cases
    )
    return (;
        id=run.id,
        material_id=run.material_id,
        source_path=spec["path"],
        source_sha256=run.provenance.sha256,
        clock_scale,
        windows,
        cases,
    )
end

is_positive_rate(diagnostic) = diagnostic.status === :positive_rate
log_rate(diagnostic) = log(diagnostic.conversion_rate_min_inv::Float64)
span(values) = maximum(values) - minimum(values)

function friedman_rate_evidence(diagnostic)
    evidence = Dict{String,Any}(
        "status" => String(diagnostic.status),
        "points_before" => diagnostic.points_before,
        "points_after" => diagnostic.points_after,
        "window_points" => diagnostic.window_points,
    )
    if diagnostic.window_source_lines !== nothing
        evidence["source_lines"] = collect(diagnostic.window_source_lines)
        evidence["time_min"] = collect(diagnostic.window_time_min)
        evidence["conversion_range"] = collect(diagnostic.window_conversion_range)
    end
    if diagnostic.conversion_rate_min_inv !== nothing
        evidence["conversion_rate_min_inv"] = diagnostic.conversion_rate_min_inv
        evidence["residual_rmse"] = diagnostic.residual_rmse
    end
    return evidence
end

function friedman_point_evidence(diagnostics)
    d = first(diagnostics)
    point = Dict{String,Any}(
        "alpha" => d.alpha,
        "crossing_status" => String(d.crossing_status),
        "rates_by_setting_order" => friedman_rate_evidence.(diagnostics),
    )
    if d.crossing_source_lines !== nothing
        point["crossing_source_lines"] = collect(d.crossing_source_lines)
        point["crossing_time_min"] = d.crossing_time_min
        point["crossing_temperature_K"] = d.crossing_temperature_K
    end
    if all(is_positive_rate, diagnostics)
        point["within_run_log_rate_span"] = span(log_rate.(diagnostics))
    end
    return point
end

function friedman_run_evidence(run)
    return Dict(
        "id" => run.id,
        "material_id" => run.material_id,
        "source_path" => run.source_path,
        "source_sha256" => run.source_sha256,
        "clock_scale_min" => run.clock_scale,
        "half_windows_min_by_setting_order" => run.windows,
        "cases" =>
            Dict(id => friedman_point_evidence.(points) for (id, points) in run.cases),
    )
end

function friedman_group_point(runs, case_id, alpha, index, setting_ids)
    diagnostics = [run.cases[case_id][index] for run in runs]
    positive_ids = Dict(
        id => [run.id for (run, d) in zip(runs, diagnostics) if is_positive_rate(d[j])] for
        (j, id) in pairs(setting_ids)
    )
    all_positive = all(ids -> length(ids) == length(runs), values(positive_ids))
    evidence = Dict{String,Any}(
        "alpha" => alpha,
        "all_settings_positive" => all_positive,
        "positive_run_ids_by_setting" => positive_ids,
    )
    if all_positive
        logs = [log_rate.(d) for d in diagnostics]
        evidence["max_within_run_log_rate_span"] = maximum(span, logs)
        evidence["max_cross_run_log_response_contrast_span"] = maximum(
            span(logs[i] .- logs[j]) for i in eachindex(logs) for j in (i + 1):length(logs)
        )
    end
    return evidence
end

function friedman_case_evidence(runs, case_id, grid, setting_ids)
    points = [
        friedman_group_point(runs, case_id, alpha, i, setting_ids) for
        (i, alpha) in pairs(grid)
    ]
    return Dict{String,Any}(
        "id" => case_id,
        "all_settings_positive_alpha" =>
            [p["alpha"] for p in points if p["all_settings_positive"]],
        "points" => points,
    )
end

function friedman_group_evidence(runs, grid, setting_ids, omitted_id)
    material = first(runs).material_id
    subset = filter(run -> run.id != omitted_id, runs)
    cases = map(sort!(collect(keys(first(runs).cases)))) do id
        evidence = friedman_case_evidence(runs, id, grid, setting_ids)
        if length(subset) < length(runs)
            omission = friedman_case_evidence(subset, id, grid, setting_ids)
            omission["run_ids"] = [r.id for r in subset]
            omission["omitted_run_id"] = omitted_id
            evidence["omit_pp_10"] = omission
        end
        return evidence
    end
    return Dict(
        "material_id" => material, "run_ids" => [r.id for r in runs], "cases" => cases
    )
end

function audit_public_tg_friedman(output=joinpath(TG_ROOT, FRIEDMAN_REPORT))
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    conversion_protocol = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    protocol = TOML.parsefile(joinpath(TG_ROOT, FRIEDMAN_PROTOCOL))
    protocol["schema_version"] == 1 || error("unsupported Friedman protocol schema")
    specs = filter(
        spec -> spec["group"] in ("experimental_tg", "synthetic_tg"), manifest["runs"]
    )
    length(specs) == 18 || error("expected the 18 pinned pICNIK TG runs")
    runs = [friedman_run(spec, manifest, conversion_protocol, protocol) for spec in specs]
    materials = sort!(unique([run.material_id for run in runs]))
    grid = conversion_protocol["alpha_grid"]
    setting_ids = [setting["id"] for setting in protocol["settings"]]
    groups = [
        friedman_group_evidence(
            filter(run -> run.material_id == material, runs),
            grid,
            setting_ids,
            protocol["omitted_polypropylene_run"],
        ) for material in materials
    ]
    evidence = Dict(
        "schema_version" => 1,
        "raw_import_manifest" => TG_IMPORT_MANIFEST,
        "conversion_protocol" => TG_PROTOCOL,
        "friedman_protocol" => FRIEDMAN_PROTOCOL,
        "comparison_case" => conversion_protocol["comparison_case"],
        "alpha_grid" => grid,
        "settings" => protocol["settings"],
        "clock_scale_alpha" => protocol["clock_scale_alpha"],
        "minimum_points_per_side" => protocol["minimum_points_per_side"],
        "groups" => groups,
        "runs" => friedman_run_evidence.(runs),
    )
    open(output, "w") do file
        return TOML.print(file, evidence; sorted=true)
    end
    return evidence
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("expected at most one report output path")
    audit_public_tg_friedman(ARGS...)
end
