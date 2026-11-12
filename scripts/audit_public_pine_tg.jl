using IsoconversionalAnalysis
using IsoconversionalAnalysis: SampledTemperatureProgram
using SHA
using Statistics
using TOML

include("audit_ictac_tg.jl")
module MassPrediction
include("audit_public_pp_prediction.jl")
end

const PINE_ROOT = normpath(joinpath(@__DIR__, ".."))
const PINE_PROTOCOL = "config/public_pine_tg.toml"
const PINE_REPORT = "docs/audits/public_pine_tg.toml"

load_pine_toml(path) = TOML.parsefile(joinpath(PINE_ROOT, path))
digest(path) = bytes2hex(open(sha256, path))

function verified_pine_file(path, expected)
    local_path = joinpath(PINE_ROOT, path)
    digest(local_path) == expected || error("Source hash mismatch for $path")
    return local_path
end

function verify_pine_sources(manifest)
    archive = manifest["archive_path"]
    archive_present = isfile(joinpath(PINE_ROOT, archive))
    if archive_present
        path = verified_pine_file(archive, manifest["archive_sha256"])
        filesize(path) == manifest["archive_bytes"] || error("Unexpected archive size")
    end
    for file in manifest["files"]
        path = verified_pine_file(file["path"], file["sha256"])
        filesize(path) == file["bytes"] || error("Unexpected file size for $path")
        haskey(file, "github_path") || continue
        snapshot = verified_pine_file(file["github_path"], file["github_sha256"])
        replace(read(path, String), "\r\n" => "\n") ==
            replace(read(snapshot, String), "\r\n" => "\n") ||
            error("TG values differ from the pinned snapshot in $path")
    end
    return archive_present
end

function pine_window(values, limits, id)
    hits = findall(x -> first(limits) <= x <= last(limits), values)
    length(hits) >= 2 || error("Insufficient samples in $limits for $id")
    rows = first(hits):last(hits)
    hits == collect(rows) || error("Discontinuous window $limits for $id")
    return rows
end

function pine_conversion(run, protocol; alternate_initial=false, alternate_terminal=false)
    T = temperature_K(run)
    t = run.time
    initial = protocol[alternate_initial ? "alternate_initial_K" : "initial_reference_K"]
    tail = protocol[alternate_terminal ? "alternate_terminal_window_s" : "terminal_window_s"]
    start = findfirst(>=(protocol["conversion_start_K"]), T)
    stop = searchsortedfirst(t, last(t) - tail)
    return tg_conversion(
        run,
        :mass;
        segment=only(run.segments).index,
        reaction_rows=start:length(t),
        initial_reference=pine_window(T, initial, run.id),
        terminal_reference=stop:length(t),
    )
end

function measured_pine_ramp(run, protocol)
    T = temperature_K(run)
    rows = pine_window(T, protocol["measured_ramp_window_K"], run.id)
    t = time_min(run)[rows]
    centered = t .- mean(t)
    return sum(centered .* (T[rows] .- mean(T[rows]))) / sum(abs2, centered)
end

function pine_prediction_input(run, protocol; held_out=false)
    T = temperature_K(run)
    t = time_min(run)
    start = findfirst(>=(protocol["prediction_start_K"]), T)
    rows = start:length(t)
    history = SampledTemperatureProgram(t[rows] .- t[start], T[rows])
    mass = signal_values(run, :mass)[rows]
    program = MassPrediction.PPProgram(history, first(mass))
    reaction = if held_out
        lo, hi = protocol["measured_hold_window_K"]
        previous = findlast(K -> !(lo <= K <= hi), history.temperature_K)
        (previous + 1):length(rows)
    else
        endpoint = findfirst(>(protocol["prediction_training_reaction_upper_K"]), T[rows])
        1:(endpoint - 1)
    end
    length(reaction) >= 2 || error("Missing scoring interval for $(run.id)")
    return MassPrediction.PPObservation(
        run.id, run.provenance.sha256, program, mass, reaction, nothing
    )
end

function pine_hold_predictions(main, holds, protocol)
    defaults = load_pine_toml(MassPrediction.PP_PROTOCOL)
    settings = merge(defaults, protocol["prediction"])
    training = [pine_prediction_input(r, protocol) for (_, r) in main]
    println("Fitting the constant-energy baseline to $(length(training)) ramps.")
    fitted = MassPrediction.fit_pp(training, nothing, settings)
    training_ids = getfield.(training, :id)
    predictions = map(holds) do (_, run)
        held = pine_prediction_input(run, protocol; held_out=true)
        record = MassPrediction.fold_record(
            "all_593K_holds", held, training, :constant_energy, fitted, Float64[], settings
        )
        predicted = pop!(record, "predicted_mass_mg")
        record["hold"] = pop!(record, "reaction_window")
        t = held.program.history.time_min
        start = first(held.reaction_rows)
        duration = last(t) - t[start]
        duration >= maximum(protocol["hold_checkpoints_min"]) || error("Truncated hold")
        last_hour = searchsortedfirst(t, last(t) - 60):length(t)
        record["hold_duration_min"] = duration
        record["recorded_hold_start_min_after_conditioning"] = t[start]
        record["conditioning_temperature_K"] = first(held.program.history.temperature_K)
        record["final_hour"] = MassPrediction.score_pp(predicted, held, last_hour)
        record["checkpoints"] = map(protocol["hold_checkpoints_min"]) do elapsed
            row = searchsortedfirst(t, t[start] + elapsed)
            Dict(
                "elapsed_hold_min" => t[row] - t[start],
                "measured_mg" => held.mass_mg[row],
                "predicted_mg" => predicted[row],
            )
        end
        return record
    end
    return Dict(
        "training_ids" => training_ids,
        "held_out_ids" => [p["held_out_run"] for p in predictions],
        "settings" => settings,
        "hold_scores" => predictions,
    )
end

function pine_mass_record(spec, run, protocol, group_key)
    curve = pine_conversion(run, protocol)
    residue = 100 * curve.terminal_reference.mass / curve.initial_reference.mass
    crossings = map(protocol["mass_series_alpha"]) do alpha
        support = tg_conversion_crossings(curve, alpha)
        record = Dict{String,Any}("alpha" => alpha, "status" => String(support.status))
        if support.status === :unique_upward
            record["temperature_K"] = only(support.crossings).temperature_K
        end
        return record
    end
    return Dict(
        "id" => run.id,
        "sha256" => spec["sha256"],
        "sample_mass_mg" => spec["sample_mass_mg"],
        group_key => spec[group_key],
        "measured_ramp_K_per_min" => measured_pine_ramp(run, protocol),
        "initial_mass_mg" => curve.initial_reference.mass,
        "terminal_mass_mg" => curve.terminal_reference.mass,
        "residue_percent" => residue,
        "crossings" => crossings,
    )
end

function pine_repeatability(entries, group_key, group_values, protocol)
    runs = [pine_mass_record(spec, run, protocol, group_key) for (spec, run) in entries]
    groups = [
        pine_group_summary(value, filter(r -> r[group_key] == value, runs)) for
        value in group_values
    ]
    return Dict("group_key" => group_key, "runs" => runs, "groups" => groups)
end

function pine_group_summary(group, runs)
    residues = [r["residue_percent"] for r in runs]
    crossing_summary = map(eachindex(first(runs)["crossings"])) do index
        crossings = [r["crossings"][index] for r in runs]
        temperatures = [r["temperature_K"] for r in crossings if haskey(r, "temperature_K")]
        record = Dict{String,Any}(
            "alpha" => first(crossings)["alpha"], "supported_repeats" => length(temperatures)
        )
        isempty(temperatures) || (record["temperature_mean_K"] = mean(temperatures))
        length(temperatures) < 2 || (record["temperature_sd_K"] = std(temperatures))
        return record
    end
    return Dict(
        "group" => group,
        "run_ids" => [r["id"] for r in runs],
        "residue_mean_percent" => mean(residues),
        "residue_sd_percent" => std(residues),
        "crossings" => crossing_summary,
    )
end

function pine_assessment(main, protocol)
    config = load_pine_toml(ICTACTGAudit.ASSESSMENT_PROTOCOL)
    config["evidence"]["pine"] = Dict("source" => "TUBS deposit README, main ramps and mass series.")
    runs = last.(main)
    rates = [s["nominal_rate_K_per_min"] for (s, _) in main]
    convert(; kwargs...) = [pine_conversion(r, protocol; kwargs...) for r in runs]
    cohort = (;
        id="pine", label="TUBS white-pine TG", synthetic=false, rates, curves=convert(),
        omit="",
        metadata=[
            Dict(
                "source_sample_mass_mg" => s["sample_mass_mg"],
                "measured_ramp_K_per_min" => measured_pine_ramp(run, protocol),
            ) for (s, run) in main
        ],
        references=Dict(
            "earlier_initial" => convert(; alternate_initial=true),
            "longer_terminal_window" => convert(; alternate_terminal=true),
        ),
    )
    evidence = [
        ICTACTGAudit.check("metadata", "2014 — sample and acquisition metadata", "met",
            "Source README reports instrument, geometry, mass, purge gas, crucible and thermal programs."),
        ICTACTGAudit.check("calibration", "2014 — calibration, blank and thermal lag", "partial",
            "README documents baseline correction and temperature/heat-flow calibration. Blank traces and independent thermal-lag measurements are unavailable."),
        ICTACTGAudit.check("repeatability", "2014 — repeatability and sample-mass sensitivity", "partial",
            "Each main rate has three repeats. The controlled mass series covers 50 K/min and cannot establish mass independence over the full range."),
        ICTACTGAudit.check("prediction", "2011 §§8–9, 2020 §2.8 — prediction on excluded programs", "met",
            "All three 593 K programs were scored together using the frozen ramp-only baseline. The status records that the test was performed."),
        ICTACTGAudit.check("components", "2020 §§2, 5 — justify a multistep model", "missing",
            "The benchmark evaluates apparent energy and a predictive baseline. Component assignments require separate evidence."),
    ]
    return ICTACTGAudit.assess_cohort(
        cohort, config, ICTACTGAudit.assessment_options(config); experiment_evidence=evidence
    )
end

function audit_public_pine_tg(output)
    protocol = load_pine_toml(PINE_PROTOCOL)
    manifest = load_pine_toml(protocol["manifest"])
    archive_present = verify_pine_sources(manifest)
    mapping = ThermalColumnMapping(
        ThermalColumnSpec("Temperature (K)", :K),
        [ThermalSignalSpec(:mass, ThermalColumnSpec("Mass (mg)", :mg), :mass)];
        time=ThermalColumnSpec("Time (s)", :s),
    )
    entries = map(filter(s -> haskey(s, "group"), manifest["files"])) do spec
        run = read_thermal_delimited(
            joinpath(PINE_ROOT, spec["path"]), mapping;
            id=spec["id"], material_id="white_pine", expected_sha256=spec["sha256"],
        )
        length(run.time) == spec["rows"] || error("Row count changed for $(run.id)")
        return spec, run
    end
    main = sort(filter(e -> first(e)["group"] == "main_ramp", entries);
        by=e -> (first(e)["nominal_rate_K_per_min"], first(e)["repeat"]))
    mass_series = sort(filter(e -> first(e)["group"] == "mass_series", entries);
        by=e -> (first(e)["nominal_mass_mg"], first(e)["repeat"]))
    holds = sort(filter(e -> first(e)["group"] == "hold", entries); by=e -> first(e)["repeat"])
    n = protocol["repeats_per_group"]
    [s["nominal_rate_K_per_min"] for (s, _) in main] ==
        repeat(protocol["main_rates_K_per_min"]; inner=n) || error("Invalid main ramp roster")
    [s["nominal_mass_mg"] for (s, _) in mass_series] ==
        repeat(protocol["mass_groups_mg"]; inner=n) || error("Invalid mass series roster")
    all(e -> first(e)["nominal_rate_K_per_min"] == protocol["mass_series_rate_K_per_min"], mass_series) ||
        error("Unexpected rate in the sample-mass series")
    length(holds) == n && all(e -> first(e)["nominal_hold_K"] == protocol["hold_temperature_K"], holds) ||
        error("Invalid holdout roster")
    println("Verified and imported $(length(entries)) TG files.")
    main_repeatability = pine_repeatability(
        main, "nominal_rate_K_per_min", protocol["main_rates_K_per_min"], protocol
    )
    mass_sensitivity = pine_repeatability(
        mass_series, "nominal_mass_mg", protocol["mass_groups_mg"], protocol
    )
    predictions = pine_hold_predictions(main, holds, protocol)
    println("Evaluating apparent energy and ICTAC evidence on the main ramps.")
    assessment = pine_assessment(main, protocol)
    inputs = [PINE_PROTOCOL, protocol["manifest"], ICTACTGAudit.ASSESSMENT_PROTOCOL,
        ICTACTGAudit.PublicTGEnergy.ENERGY_PROTOCOL, ICTACTGAudit.PublicTGEnergy.FRIEDMAN_PROTOCOL,
        MassPrediction.PP_PROTOCOL, "config/analysis_defaults.toml", "scripts/audit_public_pine_tg.jl",
        "scripts/audit_ictac_tg.jl", "scripts/audit_public_pp_prediction.jl"]
    report = Dict(
        "schema_version" => 1,
        "source_doi" => manifest["source_doi"],
        "archive_sha256" => manifest["archive_sha256"],
        "archive_rechecked_for_this_run" => archive_present,
        "github_revision" => manifest["comparison_revision"],
        "files_verified" => length(manifest["files"]),
        "tg_files_imported" => length(entries),
        "tg_rows_imported" => sum(length(run.time) for (_, run) in entries),
        "input_sha256" => Dict(path => digest(joinpath(PINE_ROOT, path)) for path in inputs),
        "assessment" => assessment,
        "main_repeatability" => main_repeatability,
        "mass_sensitivity" => mass_sensitivity,
        "prediction" => predictions,
    )
    open(output, "w") do io
        TOML.print(io, report; sorted=true)
    end
    println("Wrote $output")
    return report
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("Expected at most one output path")
    audit_public_pine_tg(isempty(ARGS) ? joinpath(PINE_ROOT, PINE_REPORT) : abspath(only(ARGS)))
end
