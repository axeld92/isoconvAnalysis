using IsoconversionalAnalysis
using TOML

module PublicTGRampInputs
include("audit_public_tg_conversion.jl")
end

using .PublicTGRampInputs:
    TG_ROOT, TG_IMPORT_MANIFEST, TG_PROTOCOL, import_tg, tg_case_curves, require_check

const RAMP_REPORT = "docs/audits/public_tg_ramp_approximations.toml"
const RAMP_PROBE_ENERGIES = [60.0, 120.0, 240.0]
const RAMP_SENSITIVITY_ALPHA = (0.10, 0.90)
const RAMP_PP_OMISSION = "picnik_PP_10_B3"

function ramp_probe_record(probe)
    return Dict{String,Any}(
        "energy_kJ_per_mol" => probe.energy_kJ_per_mol,
        "crossing_E_over_RT" => probe.crossing_E_over_RT,
        "log_measured_integral_min" => probe.log_measured_integral_min,
        "log_affine_integral_min" => probe.log_affine_integral_min,
        "measured_vs_affine_log_ratio" => probe.measured_vs_affine_log_ratio,
        "max_quadrature_log_difference" => probe.max_quadrature_log_difference,
        "omitted_start_fraction" => probe.omitted_start_fraction,
        "approximation_log_residuals" =>
            Dict(String(k) => v for (k, v) in pairs(probe.approximation_log_residuals)),
    )
end

function ramp_record(id, diagnostic)
    result = Dict{String,Any}(
        "id" => id,
        "status" => String(diagnostic.status),
        "crossing_status" => String(diagnostic.crossing.status),
        "start_row" => diagnostic.start_row,
        "start_source_line" => diagnostic.start_source_line,
        "start_temperature_K" => diagnostic.start_temperature_K,
        "start_time_min" => diagnostic.start_time_min,
    )
    if diagnostic.crossing.status === :unique_upward
        event = only(diagnostic.crossing.crossings)
        result["crossing_temperature_K"] = event.temperature_K
        result["crossing_time_min"] = event.time_min
        result["crossing_source_lines"] = collect(event.source_lines)
    end
    if diagnostic.duration_min !== nothing
        for field in (
            :duration_min,
            :secant_K_per_min,
            :OLS_K_per_min,
            :affine_max_abs_residual_K,
            :affine_RMSE_K,
            :falling_temperature_steps,
            :flat_temperature_steps,
            :cumulative_cooling_K,
            :raw_row_count,
        )
            result[String(field)] = getproperty(diagnostic, field)
        end
        result["half_secants_K_per_min"] = collect(diagnostic.half_secants_K_per_min)
        result["step_slope_range_K_per_min"] = collect(
            diagnostic.step_slope_range_K_per_min
        )
    end
    result["probes"] = ramp_probe_record.(diagnostic.probes)
    return result
end

function ramp_span(values)
    lo, hi = extrema(values)
    return hi - lo
end

function ramp_group_probe(records, i)
    probes = [record["probes"][i] for record in records]
    ratios = [p["measured_vs_affine_log_ratio"] for p in probes]
    x_range = extrema(p["crossing_E_over_RT"] for p in probes)
    shape_spans = Dict(
        method => ramp_span(p["approximation_log_residuals"][method] for p in probes) for
        method in ("kas", "fwo", "starink")
    )
    return Dict{String,Any}(
        "energy_kJ_per_mol" => RAMP_PROBE_ENERGIES[i],
        "crossing_E_over_RT_range" => collect(x_range),
        "fwo_all_crossings_within_Doyle_range" =>
            20 <= first(x_range) <= last(x_range) <= 60,
        "max_absolute_measured_vs_affine_log_ratio" => maximum(abs, ratios),
        "measured_vs_affine_log_ratio_span" => ramp_span(ratios),
        "maximum_omitted_start_fraction" =>
            maximum(p["omitted_start_fraction"] for p in probes),
        "approximation_log_residual_span" => shape_spans,
    )
end

function ramp_group_summary(records)
    supported = filter(r -> r["status"] == "evaluated", records)
    common = length(supported) == length(records)
    result = Dict{String,Any}(
        "all_prefixes_evaluated" => common,
        "evaluated_run_ids" => [record["id"] for record in supported],
    )
    if common
        result["max_affine_residual_K"] = maximum(
            r["affine_max_abs_residual_K"] for r in records
        )
        result["probes"] = [
            ramp_group_probe(records, i) for i in eachindex(RAMP_PROBE_ENERGIES)
        ]
    end
    return result
end

function ramp_public_inputs(manifest, protocol)
    return [
        begin
            run = import_tg(spec, manifest)
            selected, curves, reference = tg_case_curves(run, spec, protocol)
            require_check(
                reference.time_min !== nothing, run.id, "recorded clock required"
            )
            warm_lower_C = first(protocol["windows"][spec["group"]]["reaction_C"])
            warm_start_row = findfirst(>=(warm_lower_C), selected.temperature)
            require_check(warm_start_row !== nothing, run.id, "declared warm start absent")
            (;
                id=run.id,
                material=run.material_id,
                source_path=spec["path"],
                source_sha256=run.provenance.sha256,
                full_record_OLS_K_per_min=only(thermal_segment_diagnostics(selected)).slope_K_per_min,
                warm_start_row,
                warm_lower_C,
                curves,
            )
        end for
        spec in manifest["runs"] if spec["group"] in ("synthetic_tg", "experimental_tg")
    ]
end

function compare_ramp_starts!(warm, first_start)
    for (warm_run, first_run) in zip(warm["runs"], first_start["runs"])
        if warm_run["status"] == first_run["status"] == "evaluated"
            for (later, full) in zip(warm_run["probes"], first_run["probes"])
                fraction =
                    -expm1(
                        later["log_measured_integral_min"] -
                        full["log_measured_integral_min"],
                    )
                require_check(
                    fraction >= -1e-10,
                    warm_run["id"],
                    "shorter measured interval has larger integral",
                )
                later["fraction_of_full_measured_integral_before_warm_start"] = max(
                    0.0, fraction
                )
            end
        end
    end
    return nothing
end

function audit_ramp_point(runs, case, alpha, warm)
    records = map(runs) do run
        curve = run.curves[case]
        start_row = warm ? run.warm_start_row : first(curve.reaction_rows)
        diagnostic = tg_ramp_diagnostic(
            curve, alpha; start_row, probe_energies_kJ_per_mol=RAMP_PROBE_ENERGIES
        )
        record = ramp_record(run.id, diagnostic)
        if diagnostic.secant_K_per_min !== nothing
            record["prefix_secant_relative_to_full_record_OLS"] =
                diagnostic.secant_K_per_min / run.full_record_OLS_K_per_min
        end
        return record
    end
    point = ramp_group_summary(records)
    point["alpha"] = alpha
    point["runs"] = records
    if any(r -> r["id"] == RAMP_PP_OMISSION, records)
        point["omit_PP_10_K_per_min"] = ramp_group_summary(
            filter(r -> r["id"] != RAMP_PP_OMISSION, records)
        )
    end
    return point
end

function audit_ramp_case(runs, case, protocol)
    baseline = case == protocol["comparison_case"]
    alphas = if baseline
        protocol["alpha_grid"]
    else
        filter(
            a -> first(RAMP_SENSITIVITY_ALPHA) <= a <= last(RAMP_SENSITIVITY_ALPHA),
            protocol["alpha_grid"],
        )
    end
    first_points = [audit_ramp_point(runs, case, alpha, false) for alpha in alphas]
    starts = [Dict("id" => "selected_first_row", "points" => first_points)]
    if baseline
        warm_points = [audit_ramp_point(runs, case, alpha, true) for alpha in alphas]
        foreach(compare_ramp_starts!, warm_points, first_points)
        push!(starts, Dict("id" => "warm_start", "points" => warm_points))
    end
    return Dict("id" => case, "program_starts" => starts)
end

function audit_public_tg_ramp(output)
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    protocol = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    all_runs = ramp_public_inputs(manifest, protocol)
    require_check(length(all_runs) == 18, "pICNIK", "expected 18 public TG runs")
    groups = map(sort!(unique(run.material for run in all_runs))) do material
        runs = filter(r -> r.material == material, all_runs)
        cases = sort!(collect(keys(first(runs).curves)))
        return Dict(
            "material_id" => material,
            "run_ids" => [r.id for r in runs],
            "sources" => [
                Dict(
                    "id" => r.id,
                    "path" => r.source_path,
                    "sha256" => r.source_sha256,
                    "full_record_OLS_K_per_min" => r.full_record_OLS_K_per_min,
                ) for r in runs
            ],
            "warm_start_temperature_threshold_C" => first(runs).warm_lower_C,
            "cases" => [audit_ramp_case(runs, case, protocol) for case in cases],
        )
    end
    run_records = [
        r for g in groups for c in g["cases"] for s in c["program_starts"] for
        p in s["points"] for r in p["runs"]
    ]
    probes = [p for r in run_records for p in r["probes"]]
    max_refinement_difference = maximum(p["max_quadrature_log_difference"] for p in probes)
    require_check(
        max_refinement_difference <= 1e-8, "pICNIK", "4- and 8-node quadrature differ"
    )
    result = Dict(
        "schema_version" => 1,
        "raw_import_manifest" => TG_IMPORT_MANIFEST,
        "TG_protocol" => TG_PROTOCOL,
        "baseline_case" => protocol["comparison_case"],
        "probe_energies_kJ_per_mol" => RAMP_PROBE_ENERGIES,
        "number_of_run_prefix_checks" => length(run_records),
        "number_of_evaluated_run_prefixes" =>
            count(r -> r["status"] == "evaluated", run_records),
        "number_of_integral_probes" => length(probes),
        "max_quadrature_log_difference" => max_refinement_difference,
        "groups" => groups,
    )
    open(output, "w") do io
        return TOML.print(io, result; sorted=true)
    end
    println(
        "Wrote ",
        output,
        " with ",
        length(run_records),
        " prefix checks and ",
        length(probes),
        " integral probes.",
    )
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("Expected at most one report output path")
    output = isempty(ARGS) ? joinpath(TG_ROOT, RAMP_REPORT) : abspath(only(ARGS))
    audit_public_tg_ramp(output)
end
