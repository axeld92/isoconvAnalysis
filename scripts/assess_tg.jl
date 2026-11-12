module TGAssessment

using IsoconversionalAnalysis
using IsoconversionalAnalysis: SampledTemperatureProgram
using SHA
using Statistics
using TOML

include("audit_ictac_tg.jl")
module MassPrediction
    include("audit_public_pp_prediction.jl")
end
include("render_tg_assessment.jl")

invalid(message) = throw(ArgumentError(message))

function validate_programs(programs, specs)
    ids = [s["id"] for s in specs]
    allunique(ids) || invalid("Duplicate run identifiers")
    allunique(normpath.(getindex.(specs, "path"))) || invalid("Duplicate source paths")
    hashes = [lowercase(s["sha256"]) for s in specs if haskey(s, "sha256")]
    allunique(hashes) || invalid("Duplicate source hashes")
    allunique(getindex.(programs, "id")) || invalid("Duplicate program identifiers")

    assigned = String[]
    calibration_rates = Set{Float64}()
    holdout_rates = Set{Float64}()
    for p in programs
        id, role, kind, runs = p["id"], p["role"], p["kind"], p["runs"]
        role in ("calibration", "mass_series", "holdout") || invalid("Unknown role for $id")
        kind in ("ramp", "hold") || invalid("Unknown program kind for $id")
        kind == "hold" &&
            role != "holdout" &&
            invalid("Hold $id must be assigned to holdout")
        isempty(runs) && invalid("Program $id has no runs")
        length(runs) == get(p, "expected_repeats", length(runs)) ||
            invalid("Unexpected repeat count for $id")
        all(in(ids), runs) || invalid("Unknown run identifier in $id")
        append!(assigned, runs)
        if role == "mass_series"
            mass = get(p, "nominal_mass_mg", NaN)
            isfinite(mass) && mass > 0 || invalid("Declare a positive nominal mass for $id")
        end
        if haskey(p, "nominal_rate_K_per_min")
            rate = p["nominal_rate_K_per_min"]
            kind == "ramp" && isfinite(rate) && rate > 0 || invalid("Invalid rate for $id")
            role == "calibration" && push!(calibration_rates, rate)
            role == "holdout" && push!(holdout_rates, rate)
        elseif kind == "ramp" && role == "holdout"
            invalid("Declare the nominal rate for withheld ramp $id")
        end
    end
    allunique(assigned) || invalid("A run may belong to only one program")
    any(p -> p["role"] == "calibration", programs) || invalid("Declare calibration ramps")
    if !isempty(holdout_rates)
        all(
            p -> p["role"] != "calibration" || haskey(p, "nominal_rate_K_per_min"), programs
        ) || invalid("Declare all calibration ramp rates when withholding a ramp")
        isempty(intersect(calibration_rates, holdout_rates)) ||
            invalid("A ramp rate may not be split across calibration and holdout")
    end
    return nothing
end

function read_tg_roster(config, config_dir)
    programs = config["programs"]
    manifest = get(config, "manifest", nothing)
    if manifest === nothing
        specs = config["runs"]
    else
        haskey(config, "runs") && invalid("Use either inline runs or a manifest")
        manifest = normpath(joinpath(config_dir, manifest))
        wanted = Set(id for p in programs for id in p["runs"])
        specs = filter(
            TOML.parsefile(manifest)[get(config, "manifest_table", "runs")]
        ) do spec
            return get(spec, "id", nothing) in wanted
        end
    end
    validate_programs(programs, specs)
    specs_by_id = Dict(s["id"] => s for s in specs)
    for p in programs, id in p["runs"]
        spec = specs_by_id[id]
        haskey(spec, "nominal_hold_K") &&
            p["kind"] != "hold" &&
            invalid("Manifest identifies $id as a hold")
        if haskey(spec, "nominal_rate_K_per_min") && haskey(p, "nominal_rate_K_per_min")
            spec["nominal_rate_K_per_min"] == p["nominal_rate_K_per_min"] ||
                invalid("Program and manifest rates differ for $id")
        end
    end
    return specs_by_id, manifest
end

function import_tg_runs(config, config_dir, specs)
    columns = config["columns"]
    mapping = ThermalColumnMapping(
        ThermalColumnSpec(columns["temperature"], Symbol(columns["temperature_unit"])),
        [
            ThermalSignalSpec(
                :mass,
                ThermalColumnSpec(columns["mass"], Symbol(columns["mass_unit"])),
                :mass,
            ),
        ];
        time=ThermalColumnSpec(columns["time"], Symbol(columns["time_unit"])),
        segment_column=get(columns, "segment", nothing),
    )
    root = normpath(joinpath(config_dir, get(config, "data_root", ".")))
    runs = Dict{String,ThermalRun}()
    paths, hashes = Set{String}(), Set{String}()
    for p in config["programs"], id in p["runs"]
        spec = specs[id]
        path = realpath(joinpath(root, spec["path"]))
        path in paths && invalid("Runs share source file $path")
        push!(paths, path)
        run = read_thermal_delimited(
            path,
            mapping;
            id,
            material_id=get(config, "material_id", config["id"]),
            delimiter=only(get(columns, "delimiter", ",")),
            metadata_delimiter=only(
                get(columns, "metadata_delimiter", get(columns, "delimiter", ","))
            ),
            encoding=Symbol(get(columns, "encoding", "utf8")),
            header_prefix=get(columns, "header_prefix", nothing),
            expected_sha256=get(spec, "sha256", nothing),
        )
        length(run.time) == get(spec, "rows", length(run.time)) ||
            invalid("Row count differs for $id")
        run.provenance.sha256 in hashes &&
            invalid("Runs share raw data bytes, including $id")
        push!(hashes, run.provenance.sha256)
        segment = get(spec, "selected_segment", nothing)
        if segment === nothing
            length(run.segments) == 1 || invalid("Declare selected_segment for $id")
            segment = only(run.segments).index
        end
        runs[id] = select_thermal_segment(run, segment)
    end
    return runs
end

function tg_window(run, selector; rows=1:length(run.time))
    isempty(selector) && return rows
    length(selector) == 1 || invalid("Declare one window selector for $(run.id)")
    key = only(keys(selector))
    value = selector[key]
    T = temperature_K(run)
    if key in ("temperature_K", "temperature_C", "tail_temperature_K")
        units = key == "temperature_C" ? :C : :K
        temperatures = units == run.temperature_unit ? run.temperature : T .- 273.15
        units == :K && (temperatures = T)
        lo, hi = value
        isfinite(lo) && isfinite(hi) && lo <= hi ||
            invalid("Invalid temperature window for $(run.id)")
        if key == "tail_temperature_K"
            previous = findlast(i -> !(lo <= T[i] <= hi), rows)
            return rows[(something(previous, 0) + 1):end]
        end
        selected = filter(i -> lo <= temperatures[i] <= hi, rows)
    elseif key == "start_K"
        isfinite(value) || invalid("Invalid start temperature for $(run.id)")
        first_hit = findfirst(i -> T[i] >= value, rows)
        first_hit === nothing &&
            invalid("Start temperature is outside the record for $(run.id)")
        return rows[first_hit:end]
    elseif key == "end_K"
        isfinite(value) || invalid("Invalid end temperature for $(run.id)")
        first_above = findfirst(i -> T[i] > value, rows)
        return rows[1:(something(first_above, length(rows) + 1) - 1)]
    elseif key == "last_seconds"
        isfinite(value) && value >= 0 || invalid("Invalid tail duration for $(run.id)")
        t = run.time_unit === :s ? run.time : time_min(run) .* 60
        selected = filter(i -> t[i] >= t[last(rows)] - value, rows)
    elseif key == "rows"
        first_row, last_row = value
        first(rows) <= first_row <= last_row <= last(rows) ||
            invalid("Invalid rows for $(run.id)")
        return first_row:last_row
    else
        invalid("Unknown window selector $key")
    end
    isempty(selected) && invalid("Window $key has no observations for $(run.id)")
    length(selected) == last(selected) - first(selected) + 1 ||
        invalid("Window $key has multiple blocks for $(run.id)")
    return first(selected):last(selected)
end

function tg_curve(run, policy)
    function reference(key)
        selector = policy[key]
        if haskey(selector, "mass")
            length(selector) == 1 || invalid("Declare one $key reference for $(run.id)")
            return selector["mass"]
        end
        return tg_window(run, selector)
    end
    return tg_conversion(
        run,
        :mass;
        segment=only(run.segments).index,
        reaction_rows=tg_window(run, get(policy, "reaction", Dict())),
        initial_reference=reference("initial"),
        terminal_reference=reference("terminal"),
    )
end

function validate_output_paths(output_dir, inputs)
    # Check existing targets through symlinks and Windows path aliases before opening them.
    key(path) = Sys.iswindows() ? lowercase(realpath(path)) : realpath(path)
    protected = Set(key(path) for path in inputs)
    for name in ("tg_assessment.toml", "tg_assessment.md")
        target = joinpath(output_dir, name)
        if ispath(target) && key(target) in protected
            invalid("Report path would overwrite an input: $target")
        end
    end
    return nothing
end

function measured_rate(run, config)
    if !haskey(config, "measured_ramp_window_K")
        return only(thermal_segment_diagnostics(run)).slope_K_per_min
    end
    rows = tg_window(run, Dict("temperature_K" => config["measured_ramp_window_K"]))
    t, T = time_min(run)[rows], temperature_K(run)[rows]
    centered = t .- mean(t)
    return sum(centered .* (T .- mean(T))) / sum(abs2, centered)
end

function summarize_program(program, curves, alpha_grid)
    residue = [100 * c.terminal_reference.mass / c.initial_reference.mass for c in curves]
    crossings = map(alpha_grid) do alpha
        supports = tg_conversion_crossings.(curves, alpha)
        T = [
            only(s.crossings).temperature_K for s in supports if s.status === :unique_upward
        ]
        row = Dict{String,Any}(
            "alpha" => alpha,
            "supported_repeats" => length(T),
            "run_statuses" => String.(getproperty.(supports, :status)),
        )
        isempty(T) || (row["temperature_mean_K"] = mean(T))
        length(T) < 2 || (row["temperature_sd_K"] = std(T))
        return row
    end
    summary = Dict{String,Any}(
        "program_id" => program["id"],
        "run_ids" => program["runs"],
        "residue_mean_percent" => mean(residue),
        "crossings" => crossings,
    )
    length(residue) < 2 || (summary["residue_sd_percent"] = std(residue))
    return summary
end

function prediction_input(run, policy, window)
    start = get(policy, "start_K", nothing)
    rows =
        start === nothing ? (1:length(run.time)) : tg_window(run, Dict("start_K" => start))
    t, T = time_min(run)[rows], temperature_K(run)[rows]
    mass = signal_values(run, :mass; unit=:mg)[rows]
    first(mass) > 0 || invalid("Conditioning mass must be positive for $(run.id)")
    selected = tg_window(run, window; rows)
    length(selected) >= 2 ||
        invalid("Scoring window needs at least two points for $(run.id)")
    relative = (first(selected) - first(rows) + 1):(last(selected) - first(rows) + 1)
    program = MassPrediction.PPProgram(
        SampledTemperatureProgram(t .- first(t), T), first(mass)
    )
    return MassPrediction.PPObservation(
        run.id, run.provenance.sha256, program, mass, relative, nothing
    )
end

function frozen_prediction(config, config_dir, runs, calibration, holdouts)
    policy = config["prediction"]
    isempty(holdouts) && invalid("Prediction needs at least one complete holdout program")
    settings_path = get(policy, "settings_file", nothing)
    settings = if settings_path === nothing
        Dict()
    else
        TOML.parsefile(joinpath(config_dir, settings_path))
    end
    settings = merge(settings, get(policy, "settings", Dict()))
    training_ids = reduce(vcat, getindex.(calibration, "runs"))
    held_out_ids = reduce(vcat, getindex.(holdouts, "runs"))
    isempty(intersect(training_ids, held_out_ids)) ||
        invalid("Holdout IDs intersect calibration")
    training = [
        prediction_input(runs[id], policy, get(policy, "training_window", Dict())) for
        id in training_ids
    ]
    println("Fitting one baseline to $(length(training_ids)) calibration ramps.")
    fitted = MassPrediction.fit_pp(training, nothing, settings)
    scores = Dict{String,Any}[]
    for p in holdouts, id in p["runs"]
        held = prediction_input(runs[id], policy, get(policy, "holdout_window", Dict()))
        record = MassPrediction.fold_record(
            p["id"], held, training, :constant_energy, fitted, Float64[], settings
        )
        predicted = pop!(record, "predicted_mass_mg")
        record["evaluation"] = pop!(record, "reaction_window")
        record["program_id"] = p["id"]
        t = held.program.history.time_min
        start = t[first(held.reaction_rows)]
        record["evaluation_duration_min"] = t[last(held.reaction_rows)] - start
        if haskey(policy, "final_window_min")
            duration = policy["final_window_min"]
            0 < duration <= last(t) || invalid("Invalid final scoring window for $id")
            rows = searchsortedfirst(t, last(t) - duration):length(t)
            length(rows) >= 2 ||
                invalid("Final scoring window needs at least two points for $id")
            record["final_window"] = MassPrediction.score_pp(predicted, held, rows)
        end
        record["checkpoints"] = map(get(policy, "checkpoints_min", Float64[])) do elapsed
            0 <= elapsed <= record["evaluation_duration_min"] ||
                invalid("Checkpoint outside evaluation window for $id")
            row = searchsortedfirst(t, start + elapsed)
            return Dict(
                "elapsed_min" => t[row] - start,
                "measured_mg" => held.mass_mg[row],
                "predicted_mg" => predicted[row],
            )
        end
        push!(scores, record)
    end
    return Dict(
        "training_ids" => training_ids,
        "held_out_ids" => held_out_ids,
        "settings" => settings,
        "hold_scores" => scores,
    )
end

function experimental_checks(config, calibration, mass_series, prediction)
    evidence = get(config, "evidence", Dict())
    function check(id, recommendation; status=nothing, note="Evidence was not supplied.")
        declared = get(evidence, id, Dict())
        state = status === nothing ? get(declared, "status", "missing") : status
        state in ("met", "partial", "missing", "not met", "not applicable", "advisory") ||
            invalid("Unknown evidence status for $id")
        detail = get(declared, "note", "")
        text = status === nothing && !isempty(detail) ? detail : strip("$note $detail")
        return ICTACTGAudit.check(id, recommendation, state, text)
    end
    repeated = count(p -> length(p["runs"]) > 1, calibration)
    mass_groups = length(mass_series)
    return [
        check("metadata", "2014 — sample and acquisition metadata"),
        check("calibration", "2014 — calibration, blank and thermal lag"),
        check(
            "repeatability",
            "2014 — repeatability and sample-mass sensitivity";
            status=repeated == 0 && mass_groups == 0 ? "missing" : "partial",
            note="$repeated of $(length(calibration)) calibration programs have repeats. $mass_groups sample-mass groups were declared.",
        ),
        check(
            "prediction",
            "2011 §§8–9, 2020 §2.8 — prediction on excluded programs";
            status=prediction === nothing ? nothing : "met",
            note=if prediction === nothing
                "This config has no prediction fit."
            else
                "One frozen calibration fit was scored on $(length(prediction["held_out_ids"])) held-out runs. This records execution of the test."
            end,
        ),
        check("components", "2020 §§2, 5 — justify a multistep model"),
    ]
end

function assess_tg(config_path, output_dir)
    config_path = abspath(config_path)
    config = TOML.parsefile(config_path)
    get(config, "schema_version", 0) == 1 ||
        invalid("Unsupported TG assessment schema version")
    config_dir = dirname(config_path)
    specs, manifest = read_tg_roster(config, config_dir)
    runs = import_tg_runs(config, config_dir, specs)
    inputs = [config_path; [r.provenance.source_path for r in values(runs)]]
    manifest === nothing || push!(inputs, manifest)
    prediction_settings = get(get(config, "prediction", Dict()), "settings_file", nothing)
    prediction_settings === nothing ||
        push!(inputs, joinpath(config_dir, prediction_settings))
    validate_output_paths(output_dir, inputs)
    programs = config["programs"]
    calibration = filter(p -> p["role"] == "calibration", programs)
    mass_series = filter(p -> p["role"] == "mass_series", programs)
    holdouts = filter(p -> p["role"] == "holdout", programs)
    pairs = [(p, id) for p in calibration for id in p["runs"]]
    policy = config["conversion"]
    curves = Dict(
        id => tg_curve(runs[id], policy) for p in [calibration; mass_series] for
        id in p["runs"]
    )
    metadata = map(pairs) do (_, id)
        record = Dict{String,Any}(
            "measured_ramp_K_per_min" => measured_rate(runs[id], config)
        )
        if haskey(specs[id], "sample_mass_mg")
            record["source_sample_mass_mg"] = specs[id]["sample_mass_mg"]
        end
        return record
    end
    measured = Dict(
        id => m["measured_ramp_K_per_min"] for ((_, id), m) in zip(pairs, metadata)
    )
    program_rates = Dict{String,Float64}()
    for p in calibration
        program_rates[p["id"]] = get(p, "nominal_rate_K_per_min") do
            return mean(measured[id] for id in p["runs"])
        end
    end
    rates = [program_rates[p["id"]] for (p, _) in pairs]
    cohort = (;
        id=config["id"],
        label=config["label"],
        synthetic=false,
        rates,
        metadata,
        curves=[curves[id] for (_, id) in pairs],
        omit=get(config, "omit_run", ""),
        references=Dict(
            name => [tg_curve(runs[id], merge(policy, overrides)) for (_, id) in pairs] for
            (name, overrides) in get(policy, "variants", Dict())
        ),
    )
    ictac = ICTACTGAudit.load_toml(ICTACTGAudit.ASSESSMENT_PROTOCOL)
    ictac["evidence"][cohort.id] = Dict(
        "source" => get(get(config, "evidence", Dict()), "source", "Configured TG exports.")
    )
    prediction = if haskey(config, "prediction")
        frozen_prediction(config, config_dir, runs, calibration, holdouts)
    else
        nothing
    end
    checks = experimental_checks(config, calibration, mass_series, prediction)
    println("Assessing $(cohort.label) using $(length(pairs)) calibration runs.")
    assessment = ICTACTGAudit.assess_cohort(
        cohort, ictac, ICTACTGAudit.assessment_options(ictac); experiment_evidence=checks
    )
    alpha = get(config, "repeatability_alpha", [0.10, 0.50, 0.90])
    summaries = [
        summarize_program(p, [curves[id] for id in p["runs"]], alpha) for
        p in [calibration; mass_series]
    ]
    roster = [
        merge(Dict(k => v for (k, v) in p if k != "runs"), Dict("run_ids" => p["runs"])) for
        p in programs
    ]
    source_records = map([id for p in programs for id in p["runs"]]) do id
        run, spec = runs[id], specs[id]
        record = Dict{String,Any}(
            "id" => id,
            "path" => spec["path"],
            "sha256" => run.provenance.sha256,
            "rows" => length(run.time),
            "segment" => only(run.segments).index,
        )
        if haskey(spec, "sample_mass_mg")
            record["sample_mass_mg"] = spec["sample_mass_mg"]
        end
        return record
    end
    numerical_configs = [
        ICTACTGAudit.ASSESSMENT_PROTOCOL,
        ICTACTGAudit.PublicTGEnergy.ENERGY_PROTOCOL,
        ICTACTGAudit.PublicTGEnergy.FRIEDMAN_PROTOCOL,
        "config/analysis_defaults.toml",
    ]
    report = Dict{String,Any}(
        "schema_version" => 1,
        "config_sha256" => bytes2hex(open(sha256, config_path)),
        "numerical_config_sha256" => Dict(
            p => bytes2hex(open(sha256, joinpath(ICTACTGAudit.ROOT, p))) for
            p in numerical_configs
        ),
        "references" => ictac["references"],
        "assessment" => assessment,
        "programs" => roster,
        "repeatability" => summaries,
        "runs" => source_records,
    )
    manifest === nothing || (report["manifest_sha256"] = bytes2hex(open(sha256, manifest)))
    prediction === nothing || (report["prediction"] = prediction)
    mkpath(output_dir)
    open(joinpath(output_dir, "tg_assessment.toml"), "w") do io
        return TOML.print(io, report; sorted=true)
    end
    open(joinpath(output_dir, "tg_assessment.md"), "w") do io
        return render_tg_assessment(io, report, "tg_assessment.toml")
    end
    println("Wrote TG assessment to $(abspath(output_dir)).")
    return report
end

end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 ||
        error("Usage: julia --project=. scripts/assess_tg.jl CONFIG OUTPUT_DIR")
    TGAssessment.assess_tg(ARGS...)
end
