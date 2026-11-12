module ICTACTGAudit

using IsoconversionalAnalysis
using Printf
using SHA
using Statistics
using TOML

module PublicTGEnergy
    include("audit_public_tg_energy.jl")
end

using .PublicTGEnergy.PublicTGEnergyInputs.PublicTGFriedmanInputs: window_rows

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ASSESSMENT_PROTOCOL = "config/ictac_tg.toml"
const METHODS = (:friedman, :advanced_vyazovkin)

load_toml(path) = TOML.parsefile(joinpath(ROOT, path))
number(value) = @sprintf("%.3f", value)
function check(id, recommendation, status, evidence)
    return Dict(
        "id" => id,
        "recommendation" => recommendation,
        "status" => status,
        "evidence" => evidence,
    )
end

function assessment_options(config)
    grid = config["alpha_grid"]
    length(grid) >= 2 &&
    first(grid) == 0.05 &&
    last(grid) == 0.95 &&
    all(d -> 0 < d <= 0.05 + 8eps(Float64), diff(grid)) ||
        error("The ICTAC reporting grid must cover 0.05–0.95 in steps of at most 0.05")
    energy = load_toml(PublicTGEnergy.ENERGY_PROTOCOL)
    derivative = load_toml(PublicTGEnergy.FRIEDMAN_PROTOCOL)
    defaults = load_config(joinpath(ROOT, "config", "analysis_defaults.toml"))
    energy["alpha_grid"] = grid
    return (;
        protocol=energy,
        advanced=defaults.isoconversional.advanced_vyazovkin,
        clock_bounds=derivative["clock_scale_alpha"],
        minimum_points=derivative["minimum_points_per_side"],
        primary_derivative=only(
            filter(s -> s["id"] == energy["primary_derivative"], derivative["settings"])
        ),
        derivative_settings=derivative["settings"],
        omitted_run=derivative["omitted_polypropylene_run"],
    )
end

function program_checks(rates)
    !isempty(rates) && all(r -> isfinite(r) && r > 0, rates) ||
        throw(ArgumentError("Declare a positive rate for each program"))
    distinct = sort!(unique(rates))
    ratio = last(distinct) / first(distinct)
    detail = "$(length(distinct)) distinct heating programs across $(length(rates)) runs."
    return [
        check(
            "programs",
            "2011 §§2.3, 3 and 9 — at least three different programs",
            length(distinct) >= 3 ? "met" : "not met",
            detail,
        ),
        check(
            "spread",
            "2020 §1 — aim for roughly a decade of heating rates",
            ratio >= 10 ? "met" : "advisory",
            "Fastest/slowest declared or measured rate is $(number(ratio))×.",
        ),
    ]
end

function fit_profiles(
    curves, options, case; derivative=options.primary_derivative, methods=METHODS
)
    # Use the same derivative clock convention as the public PP workflow. Retain
    # a documented duration fallback so a missing anchor does not erase other points.
    fallback = String[]
    clocks = map(curves) do curve
        crossings = tg_conversion_crossings.(Ref(curve), options.clock_bounds)
        if all(c -> c.status === :unique_upward, crossings)
            return PublicTGEnergy.friedman_clock_scale(curve, options.clock_bounds)
        end
        push!(fallback, curve.run_id)
        return last(curve.time_min) - first(curve.time_min)
    end
    return map(collect(methods)) do method
        profile = PublicTGEnergy.fit_energy(
            method, curves, options; case, derivative, clocks
        )
        method === :friedman && (profile["settings"]["duration_fallback_runs"] = fallback)
        return profile
    end
end

function estimated(profile, bounds)
    return Dict(
        p["alpha"] => p["energy_kJ_per_mol"] for p in profile["points"] if
        p["status"] == "estimated" && first(bounds) <= p["alpha"] <= last(bounds)
    )
end

function central_range(profile, bounds)
    values_at_alpha = estimated(profile, bounds)
    values_in_order = [values_at_alpha[a] for a in sort!(collect(keys(values_at_alpha)))]
    result = Dict{String,Any}("estimated_points" => length(values_in_order))
    if length(values_in_order) >= 2 && mean(values_in_order) > 0
        result["range_over_mean_percent"] =
            100 * (maximum(values_in_order) - minimum(values_in_order)) /
            mean(values_in_order)
    end
    return result
end

function profile_shift(baseline, variant, bounds)
    a, b = estimated(baseline, bounds), estimated(variant, bounds)
    common = sort!(collect(intersect(keys(a), keys(b))))
    result = Dict{String,Any}(
        "case" => variant["case"], "method" => variant["method"], "common_alpha" => common
    )
    if !isempty(common)
        result["max_abs_shift_kJ_per_mol"] = maximum(abs(a[x] - b[x]) for x in common)
    end
    return result
end

function sensitivity_profiles(cohort, baseline, options, config)
    variants = Dict{String,Any}[]
    for id in config["derivative_sensitivity"]
        settings = only(filter(s -> s["id"] == id, options.derivative_settings))
        append!(
            variants,
            fit_profiles(
                cohort.curves, options, id; derivative=settings, methods=(:friedman,)
            ),
        )
    end
    for name in sort!(collect(keys(cohort.references)))
        append!(variants, fit_profiles(cohort.references[name], options, name))
    end
    if !isempty(cohort.omit)
        kept = filter(c -> c.run_id != cohort.omit, cohort.curves)
        length(kept) == length(cohort.curves) - 1 ||
            error("Unexpected omission for $(cohort.id)")
        append!(variants, fit_profiles(kept, options, "omit_" * cohort.omit))
    end
    summaries = map(variants) do variant
        reference = only(filter(p -> p["method"] == variant["method"], baseline))
        return profile_shift(reference, variant, config["central_alpha"])
    end
    return variants, summaries
end

function source_record(curve, rate, config; synthetic)
    path = if synthetic
        joinpath("generated", basename(curve.source_path))
    else
        relpath(curve.source_path, ROOT)
    end
    ramp = tg_ramp_diagnostic(
        curve,
        config["ramp_probe_alpha"];
        start_row=first(curve.reaction_rows),
        probe_energies_kJ_per_mol=[config["ramp_probe_energy_kJ_per_mol"]],
    )
    reference(r) = Dict(
        "value" => r.mass,
        "method" => String(r.method),
        "source_lines" => r.source_lines,
    )
    record = Dict{String,Any}(
        "run_id" => curve.run_id,
        "program_rate_K_per_min" => rate,
        "path" => replace(path, '\\' => '/'),
        "sha256" => curve.source_sha256,
        "segment" => curve.segment_index,
        "reaction_source_lines" => [first(curve.source_lines), last(curve.source_lines)],
        "mass_unit" => String(curve.mass_unit),
        "initial" => reference(curve.initial_reference),
        "terminal" => reference(curve.terminal_reference),
        "negative_mass_rows" => curve.negative_mass_count,
        "decreasing_conversion_steps" => curve.decreasing_steps,
        "outside_unit_interval_rows" => curve.outside_unit_interval_count,
        "ramp_probe_status" => String(ramp.status),
    )
    if !isempty(ramp.probes)
        record["max_affine_temperature_residual_K"] = ramp.affine_max_abs_residual_K
        record["measured_vs_affine_integral_percent"] =
            100 * expm1(only(ramp.probes).measured_vs_affine_log_ratio)
    end
    return record
end

function experiment_checks(cohort, config)
    if cohort.synthetic
        return [
            check(
                "metadata",
                "2014 — sample and acquisition metadata",
                "not applicable",
                "Constructed temperature, clock and mass observations with a known law.",
            ),
            check(
                "calibration",
                "2014 — calibration, blank and thermal lag",
                "not applicable",
                "The generator models an ideal instrument.",
            ),
            check(
                "repeatability",
                "2014 — repeatability and sample-mass sensitivity",
                "not applicable",
                "Deterministic numerical controls.",
            ),
            check(
                "prediction",
                "2011 §§8–9, 2020 §2.8 — prediction on excluded programs",
                "not applicable",
                "The control tests energy recovery.",
            ),
            check(
                "components",
                "2020 §§2, 5 — justify a multistep model",
                "not applicable",
                "The prescribed control law supplies the kinetic truth.",
            ),
        ]
    end
    evidence = config["evidence"][cohort.id]
    tires = cohort.id == "tires"
    return [
        check(
            "metadata",
            "2014 — sample and acquisition metadata",
            "partial",
            evidence["metadata"],
        ),
        check(
            "calibration",
            "2014 — calibration, blank and thermal lag",
            tires ? "partial" : "missing",
            evidence["calibration"],
        ),
        check(
            "repeatability",
            "2014 — repeatability and sample-mass sensitivity",
            "missing",
            evidence["repeatability"],
        ),
        check(
            "prediction",
            "2011 §§8–9, 2020 §2.8 — prediction on excluded programs",
            tires ? "missing" : "partial",
            evidence["prediction"],
        ),
        check(
            "components",
            "2020 §§2, 5 — justify a multistep model",
            tires ? "not met" : "missing",
            cohort.component_evidence,
        ),
    ]
end

function assess_cohort(
    cohort,
    config,
    options;
    truth=nothing,
    experiment_evidence=experiment_checks(cohort, config),
)
    checks = program_checks(cohort.rates)
    length(unique(cohort.rates)) >= 3 || error("$(cohort.id) needs three distinct programs")
    baseline = fit_profiles(cohort.curves, options, "reference")
    runs = [
        source_record(curve, rate, config; synthetic=cohort.synthetic) for
        (curve, rate) in zip(cohort.curves, cohort.rates)
    ]
    if hasproperty(cohort, :metadata)
        for (run, metadata) in zip(runs, cohort.metadata)
            run["acquisition_metadata"] = metadata
        end
    end
    variants, shifts = sensitivity_profiles(cohort, baseline, options, config)
    coverage = [
        count(p -> p["status"] == "estimated", profile["points"]) for profile in baseline
    ]
    targets = length(config["alpha_grid"])
    push!(
        checks,
        check(
            "reporting_grid",
            "2011 §3 — 0.05–0.95 with step at most 0.05",
            all(==(targets), coverage) ? "met" : "partial",
            "Friedman $(coverage[1])/$targets and advanced Vyazovkin $(coverage[2])/$targets estimated. See point statuses.",
        ),
    )
    negative = [
        r["run_id"] for
        r in runs if r["negative_mass_rows"] > 0 || r["terminal"]["value"] < 0
    ]
    push!(
        checks,
        check(
            "conversion",
            "2011 §§1, 2 and 2023 §6 — document conversion and measured residue",
            isempty(negative) ? "met" : "partial",
            if isempty(negative)
                "Run intervals and mass references are recorded."
            else
                "Intervals and references are recorded. Negative apparent mass in " *
                join(negative, ", ") *
                " limits physical interpretation of residue."
            end,
        ),
    )
    evaluated = count(r -> r["ramp_probe_status"] == "evaluated", runs)
    push!(
        checks,
        check(
            "methods",
            "2020 §4 — use differential or flexible integral estimates",
            evaluated == length(runs) ? "met" : "partial",
            "Primary fits use Friedman and recorded-time advanced Vyazovkin. The separate ramp probe covers $evaluated/$(length(runs)) runs at the declared conversion and energy.",
        ),
    )
    compared = count(s -> !isempty(s["common_alpha"]), shifts)
    push!(
        checks,
        check(
            "sensitivity",
            "2011 §2 and 2020 §1 — check stability of the reported energy",
            compared == length(shifts) ? "met" : "partial",
            "$compared/$(length(shifts)) requested comparisons have common estimated central points. These are scenario shifts.",
        ),
    )
    push!(
        checks,
        check(
            "uncertainty",
            "2020 §1 — judge energy variation against uncertainty",
            cohort.synthetic ? "not applicable" : "missing",
            if cohort.synthetic
                "Recovery is scored directly against prescribed truth."
            else
                "The central range/mean and scenario shifts are descriptive. Experimental and whole-pipeline uncertainty remain unquantified."
            end,
        ),
    )
    append!(checks, experiment_evidence)
    if truth !== nothing
        foreach(p -> PublicTGEnergy.score_control!(p, truth, options.protocol), baseline)
        push!(
            checks,
            check(
                "known_truth",
                "Numerical benchmark — recover the prescribed energy",
                "met",
                "Both primary methods passed the declared recovery tolerances.",
            ),
        )
    end
    return Dict{String,Any}(
        "id" => cohort.id,
        "label" => cohort.label,
        "checks" => checks,
        "runs" => runs,
        "baseline" => baseline,
        "sensitivity_profiles" => variants,
        "sensitivity" => shifts,
        "central_range" => Dict(
            p["method"] => central_range(p, config["central_alpha"]) for p in baseline
        ),
        "study_evidence" => if cohort.synthetic
            "Independent RK4 generator and declared truth."
        else
            config["evidence"][cohort.id]["source"]
        end,
    )
end

function polypropylene_cohort(config, options)
    protocol = load_toml(PublicTGEnergy.TG_PROTOCOL)
    runs = PublicTGEnergy.pp_runs(protocol, options)
    curves_for(key) = TGConversionResult[r.cases[key] for r in runs]
    return (;
        id="polypropylene",
        label="Polypropylene",
        synthetic=false,
        curves=curves_for(protocol["comparison_case"]),
        rates=getfield.(runs, :rate),
        omit=options.omitted_run,
        references=Dict(
            name => curves_for(name) for name in config["pp_reference_sensitivity"]
        ),
        component_evidence="No component triplet has been established for these files.",
    )
end

function tire_cohort(config)
    catalog = load_toml("config/datasets.toml")
    source = only(filter(s -> s["id"] == "dynamic_2021", catalog["sources"]))
    specs = sort!(
        filter(r -> r["waste_tire_dry_mass_percent"] == 100.0, source["experiments"]);
        by=r -> r["nominal_heating_rate_K_per_min"],
    )
    mapping = ThermalColumnMapping(
        ThermalColumnSpec("Temp./°C", :C),
        [ThermalSignalSpec(:mass, ThermalColumnSpec("Mass/%", :percent), :mass)];
        time=ThermalColumnSpec("Time/min", :min),
        segment_column="Segment",
    )
    policy = config["tire"]
    cases = map(specs) do spec
        path = normpath(joinpath(ROOT, "config", spec["raw_export_file"]))
        run = read_thermal_delimited(
            path,
            mapping;
            id=spec["id"],
            material_id="waste_tires",
            encoding=:cp1252,
            delimiter=(';'),
            header_prefix="##",
            metadata_delimiter=(':'),
        )
        ramp = only(filter(segment -> segment.label == "2", run.segments))
        selected = select_thermal_segment(run, ramp.index)
        temperature = selected.temperature
        measured = only(thermal_segment_diagnostics(selected)).slope_K_per_min
        metadata = Dict{String,Any}(
            "measured_ramp_K_per_min" => measured,
            "sample_mass_mg" => parse(Float64, run.metadata["SAMPLE MASS /mg"]),
            "instrument" => run.metadata["INSTRUMENT"],
            "calibration_file" => run.metadata["TEMPCAL"],
            "correction_file" => run.metadata["CORR. FILE"],
        )
        metadata["sample_mass_mg"] == spec["initial_mass_mg"] ||
            error("Export mass differs from catalog for $(run.id)")
        rows = window_rows(temperature, policy["reaction_C"])
        function convert(initial, terminal)
            return tg_conversion(
                run,
                :mass;
                segment=ramp.index,
                reaction_rows=rows,
                initial_reference=window_rows(temperature, policy[initial]),
                terminal_reference=window_rows(temperature, policy[terminal]),
            )
        end
        return (
            reference=convert("initial_C", "terminal_C"),
            alternative_initial=convert("alternative_initial_C", "terminal_C"),
            alternative_terminal=convert("initial_C", "alternative_terminal_C"),
            metadata=metadata,
        )
    end
    return (;
        id="tires",
        label="Current waste-tire ramps",
        synthetic=false,
        curves=getfield.(cases, :reference),
        rates=[s["nominal_heating_rate_K_per_min"] for s in specs],
        metadata=getfield.(cases, :metadata),
        omit=policy["omitted_run"],
        references=Dict(
            String(name) => getfield.(cases, name) for
            name in (:alternative_initial, :alternative_terminal)
        ),
        component_evidence="The existing component audit found inconsistent correspondence in both three- and four-peak models. It retained warnings and did not score frozen component predictions.",
    )
end

function assess_controls(config, options)
    return mktempdir() do directory
        runs = PublicTGEnergy.generate_tg_energy_controls(directory, options.protocol)
        return map(options.protocol["controls"]) do truth
            id = truth["id"]
            specs = filter(s -> s["material_id"] == id, runs)
            curves = [
                PublicTGEnergy.load_energy_control(directory, s, options.protocol) for
                s in specs
            ]
            cohort = (;
                id,
                label="Known truth — " * replace(id, '_' => ' '),
                synthetic=true,
                curves,
                rates=getindex.(specs, "rate_K_per_min"),
                omit="",
                references=Dict{String,Vector{TGConversionResult}}(),
            )
            return assess_cohort(cohort, config, options; truth)
        end
    end
end

function render_report(
    io,
    report;
    generator="scripts/audit_ictac_tg.jl",
    evidence_link="audits/ictac_tg_assessment.toml",
)
    println(io, "# ICTAC-aligned assessment of TG kinetics\n")
    println(
        io,
        "Generated by \`$generator\`. This assessment covers the ",
        "declared experiments and apparent activation energy. The accompanying ",
        "[TOML]($evidence_link) retains the numerical settings, ",
        "source hashes, conversion windows, and per-point diagnostics.\n",
    )
    println(
        io,
        "This assessment follows the ICTAC recommendations on ",
        "[2011 computation]($(report["references"]["computation"])), ",
        "[2014 data collection]($(report["references"]["data_collection"])), ",
        "[2020 multistep kinetics]($(report["references"]["multistep"])), and ",
        "[2023 decomposition]($(report["references"]["decomposition"])).\n",
    )
    println(
        io,
        "The approximately tenfold ramp spread is advisory. A reported central ",
        "energy range describes the estimates over α = 0.10–0.90. The 2020 paper's ",
        "10–20% guide depends on uncertainty. Scenario shifts are sensitivity ",
        "evidence and the fixed ramp probe is a diagnostic at one declared conversion. ",
        "Each status applies to the individual check.\n",
    )
    for cohort in report["cohorts"]
        println(io, "## $(cohort["label"])\n")
        println(io, "| Recommendation or check | Status | Evidence |")
        println(io, "|---|---|---|")
        for row in cohort["checks"]
            println(
                io,
                "| $(row["recommendation"]) | **$(row["status"])** | $(row["evidence"]) |",
            )
        end
        println(io, "\nEvidence used for this cohort. $(cohort["study_evidence"])\n")
        println(io, "| α | Friedman, kJ/mol | Advanced Vyazovkin, kJ/mol |")
        println(io, "|---:|---:|---:|")
        for (f, v) in zip((p["points"] for p in cohort["baseline"])...)
            cell(p) =
                if p["status"] == "estimated"
                    number(p["energy_kJ_per_mol"])
                else
                    replace(p["status"], '_' => ' ')
                end
            println(io, "| $(number(f["alpha"])) | $(cell(f)) | $(cell(v)) |")
        end
        println(io, "\n| Method | Central points | Central range/mean (%) |")
        println(io, "|---|---:|---:|")
        for method in sort!(collect(keys(cohort["central_range"])))
            values = cohort["central_range"][method]
            range = if haskey(values, "range_over_mean_percent")
                number(values["range_over_mean_percent"])
            else
                "—"
            end
            println(io, "| $method | $(values["estimated_points"]) | $range |")
        end
        println(
            io,
            "\n| Sensitivity case | Method | Common central points | Max shift (kJ/mol) |",
        )
        println(io, "|---|---|---:|---:|")
        for row in cohort["sensitivity"]
            shift = if haskey(row, "max_abs_shift_kJ_per_mol")
                number(row["max_abs_shift_kJ_per_mol"])
            else
                "—"
            end
            println(
                io,
                "| $(row["case"]) | $(row["method"]) | $(length(row["common_alpha"])) | $shift |",
            )
        end
        println(io)
    end
end

function main(args=ARGS)
    length(args) <= 1 ||
        error("Usage: julia --project=. scripts/audit_ictac_tg.jl [output_directory]")
    config = load_toml(ASSESSMENT_PROTOCOL)
    options = assessment_options(config)
    cohorts = assess_controls(config, options)
    for cohort in (polypropylene_cohort(config, options), tire_cohort(config))
        println("Assessing $(cohort.label).")
        push!(cohorts, assess_cohort(cohort, config, options))
    end
    inputs = [
        ASSESSMENT_PROTOCOL,
        PublicTGEnergy.ENERGY_PROTOCOL,
        PublicTGEnergy.FRIEDMAN_PROTOCOL,
        PublicTGEnergy.TG_PROTOCOL,
        PublicTGEnergy.TG_IMPORT_MANIFEST,
        "config/datasets.toml",
        "config/analysis_defaults.toml",
        "scripts/audit_ictac_tg.jl",
        "scripts/generate_tg_energy_controls.jl",
    ]
    report = Dict(
        "schema_version" => 1,
        "protocol" => ASSESSMENT_PROTOCOL,
        "references" => config["references"],
        "alpha_grid" => config["alpha_grid"],
        "central_alpha" => config["central_alpha"],
        "ramp_probe_alpha" => config["ramp_probe_alpha"],
        "ramp_probe_energy_kJ_per_mol" => config["ramp_probe_energy_kJ_per_mol"],
        "input_sha256" =>
            Dict(p => bytes2hex(open(sha256, joinpath(ROOT, p))) for p in inputs),
        "cohorts" => cohorts,
    )
    directory = isempty(args) ? joinpath(ROOT, "docs") : abspath(only(args))
    evidence = joinpath(directory, "audits", "ictac_tg_assessment.toml")
    markdown = joinpath(directory, "ictac_tg_assessment.md")
    mkpath(dirname(evidence))
    open(evidence, "w") do io
        return TOML.print(io, report; sorted=true)
    end
    open(markdown, "w") do io
        return render_report(io, report)
    end
    println("Wrote $markdown and $evidence")
    return report
end

end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    ICTACTGAudit.main()
end
