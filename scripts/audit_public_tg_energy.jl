using IsoconversionalAnalysis
using SHA
using Statistics
using TOML

include("generate_tg_energy_controls.jl")
using .TGEnergyControls: generate_tg_energy_controls

module PublicTGEnergyInputs
include("audit_public_tg_friedman.jl")
end
using .PublicTGEnergyInputs:
    TG_ROOT,
    TG_IMPORT_MANIFEST,
    TG_PROTOCOL,
    FRIEDMAN_PROTOCOL,
    import_tg,
    tg_case_curves,
    friedman_clock_scale

const ENERGY_PROTOCOL = "config/public_tg_energy.toml"
const ENERGY_REPORT = "docs/audits/public_tg_energy.toml"
const RAMP_METHODS = (:kas, :fwo, :starink)

function energy_search(base, delta)
    return AdvancedVyazovkinConfig(
        delta,
        base.integration_points,
        base.energy_bounds_kJ_per_mol,
        base.optimizer_grid_points,
        base.optimizer_tolerance_kJ_per_mol,
        base.optimizer_max_iterations,
        base.boundary_tolerance_kJ_per_mol,
        base.flat_relative_tolerance,
    )
end

function energy_point_record(point)
    s = point.support
    row = Dict{String,Any}(
        "alpha" => s.alpha,
        "status" => String(point.status),
        "run_statuses" => String.(point.run_statuses),
    )
    point.energy_kJ_per_mol === nothing ||
        (row["energy_kJ_per_mol"] = point.energy_kJ_per_mol)
    if s.common_upward
        row["crossing_temperature_K"] = [
            only(d.crossings).temperature_K for d in s.crossings
        ]
        row["crossing_source_lines"] = [
            collect(only(d.crossings).source_lines) for d in s.crossings
        ]
        row["temperature_span_K"] = s.temperature_span_K
        row["reciprocal_temperature_sxx_K_inv2"] = s.reciprocal_temperature_sxx_K_inv2
    end
    if !isempty(point.method_inputs)
        row["method_inputs"] = point.method_inputs
    end
    if point.regression !== nothing
        r = point.regression
        row["slope_K"] = r.slope
        row["log_response_residuals"] = r.residuals
        isfinite(r.r_squared) && (row["r_squared"] = r.r_squared)
    end
    if point.objective !== nothing
        n = length(point.run_statuses)
        row["objective_minimum"] = point.objective
        row["objective_excess_over_equal_integrals"] = point.objective - n * (n - 1)
        row["objective_relative_range"] = point.objective_relative_range
        row["objective_local_minima"] = point.objective_local_minima
    end
    return row
end

function fit_energy(
    method,
    curves,
    options;
    case,
    derivative=nothing,
    clocks=nothing,
    delta=options.protocol["primary_delta_alpha"],
    subdivisions=options.protocol["time_subdivisions"],
    ramp_rates=nothing,
)
    targets = options.protocol["alpha_grid"]
    settings = Dict{String,Any}()
    if method === :friedman
        half_windows = clocks .* derivative["half_window_fraction"]
        settings["derivative"] = derivative["id"]
        settings["half_windows_min"] = half_windows
        settings["degree"] = derivative["degree"]
        fit = tg_activation_energy(
            method,
            curves,
            targets;
            friedman_half_windows_min=half_windows,
            friedman_degree=derivative["degree"],
            minimum_points_per_side=options.minimum_points,
        )
    elseif method === :advanced_vyazovkin
        advanced = energy_search(options.advanced, delta)
        settings["delta_alpha"] = delta
        settings["time_subdivisions"] = subdivisions
        settings["energy_bounds_kJ_per_mol"] = collect(advanced.energy_bounds_kJ_per_mol)
        fit = tg_activation_energy(
            method, curves, targets; advanced, time_subdivisions=subdivisions
        )
    else
        settings["declared_ramp_rates_K_per_min"] = ramp_rates
        fit = tg_activation_energy(method, curves, targets; ramp_rates_K_per_min=ramp_rates)
    end
    return Dict{String,Any}(
        "case" => case,
        "method" => String(fit.method),
        "run_ids" => fit.run_ids,
        "settings" => settings,
        "points" => energy_point_record.(fit.points),
    )
end

energies(profile) = Float64[p["energy_kJ_per_mol"] for p in profile["points"]]

function require_estimated(profile)
    for point in profile["points"]
        point["status"] == "estimated" || error(
            "$(profile["case"]) $(profile["method"]) alpha=$(point["alpha"]): $(point["status"])",
        )
    end
    return profile
end

function load_energy_control(directory, spec, protocol)
    mapping = ThermalColumnMapping(
        ThermalColumnSpec("temperature_K", :K),
        [ThermalSignalSpec(:mass, ThermalColumnSpec("mass_mg", :mg), :mass)];
        time=ThermalColumnSpec("time_min", :min),
    )
    run = read_thermal_delimited(
        joinpath(directory, spec["filename"]),
        mapping;
        id=spec["id"],
        material_id=spec["material_id"],
        expected_sha256=spec["sha256"],
    )
    length(run.temperature) == protocol["control_rows"] ||
        error("control row count changed")
    return tg_conversion(
        run,
        :mass;
        segment=1,
        reaction_rows=1:length(run.temperature),
        initial_reference=protocol["control_mass_mg"][1],
        terminal_reference=protocol["control_mass_mg"][2],
    )
end

function score_control!(profile, control, protocol)
    require_estimated(profile)
    point_alpha = [p["alpha"] for p in profile["points"]]
    intercept = control["energy_intercept_kJ_per_mol"]
    slope = control["energy_slope_kJ_per_mol"]
    expected = intercept .+ slope .* point_alpha
    error_at_target = maximum(abs, energies(profile) .- expected)
    profile["truth_at_target_kJ_per_mol"] = expected
    profile["max_abs_error_at_target_kJ_per_mol"] = error_at_target
    method = profile["method"]
    limits = protocol["recovery_limit_kJ_per_mol"]
    if slope != 0 && method == "advanced_vyazovkin"
        delta = profile["settings"]["delta_alpha"]
        midpoint = expected .- slope * delta / 2
        midpoint_error = maximum(abs, energies(profile) .- midpoint)
        profile["truth_at_interval_midpoint_kJ_per_mol"] = midpoint
        profile["max_abs_difference_from_midpoint_kJ_per_mol"] = midpoint_error
        midpoint_error <= limits["variable_increment_midpoint"] ||
            error("variable-energy increment recovery failed")
    else
        limit_key = slope == 0 ? "constant_" * method : "variable_friedman"
        error_at_target <= limits[limit_key] || error(
            "$(control["id"]) $method recovery error $error_at_target exceeds $(limits[limit_key])",
        )
    end
    return profile
end

function audit_controls(options)
    return mktempdir() do directory
        manifest = generate_tg_energy_controls(directory, options.protocol)
        results = map(options.protocol["controls"]) do control
            case = control["id"]
            specs = filter(spec -> spec["material_id"] == case, manifest)
            getindex.(specs, "rate_K_per_min") ==
            options.protocol["control_rates_K_per_min"] ||
                error("control program roster changed")
            curves = [
                load_energy_control(directory, spec, options.protocol) for spec in specs
            ]
            clocks = [friedman_clock_scale(curve, options.clock_bounds) for curve in curves]
            profiles = [
                fit_energy(
                    :friedman,
                    curves,
                    options;
                    case,
                    derivative=options.primary_derivative,
                    clocks,
                ),
            ]
            increments = if control["energy_slope_kJ_per_mol"] == 0
                [options.protocol["primary_delta_alpha"]]
            else
                options.increments
            end
            for delta in increments
                push!(
                    profiles, fit_energy(:advanced_vyazovkin, curves, options; case, delta)
                )
            end
            if case == "constant_affine"
                for method in RAMP_METHODS
                    push!(
                        profiles,
                        fit_energy(
                            method,
                            curves,
                            options;
                            case,
                            ramp_rates=options.protocol["control_rates_K_per_min"],
                        ),
                    )
                end
            end
            for profile in profiles
                score_control!(profile, control, options.protocol)
            end
            if length(increments) > 1
                errors = [
                    profile["max_abs_error_at_target_kJ_per_mol"] for
                    profile in profiles[2:end]
                ]
                all(>(0), diff(errors)) || error(
                    "narrowing the increment did not improve pointwise variable-energy recovery",
                )
            end
            println("Known-truth recovery passed for $case.")
            return Dict("definition" => control, "runs" => specs, "profiles" => profiles)
        end
        return results
    end
end

function pp_runs(conversion_protocol, options)
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    runs = map(filter(spec -> spec["group"] == "experimental_tg", manifest["runs"])) do spec
        run = import_tg(spec, manifest)
        selected, cases, reference = tg_case_curves(run, spec, conversion_protocol)
        clock = friedman_clock_scale(reference, options.clock_bounds)
        rate = only(thermal_segment_diagnostics(selected)).slope_K_per_min
        return (;
            id=run.id,
            path=spec["path"],
            sha256=run.provenance.sha256,
            clock,
            rate=rate::Float64,
            cases,
            terminal_mg=reference.terminal_reference.mass,
        )
    end
    length(runs) == 5 && count(r -> r.id == options.omitted_run, runs) == 1 ||
        error("unexpected polypropylene cohort")
    return sort!(runs; by=r -> r.rate)
end

function audit_pp_cohort(cohort, runs, comparison_case, options)
    clocks = getfield.(runs, :clock)
    profiles = Dict{String,Any}[]
    for case in sort!(collect(keys(first(runs).cases)))
        curves = TGConversionResult[r.cases[case] for r in runs]
        push!(
            profiles,
            fit_energy(
                :friedman,
                curves,
                options;
                case,
                derivative=options.primary_derivative,
                clocks,
            ),
        )
        push!(profiles, fit_energy(:advanced_vyazovkin, curves, options; case))
    end
    case = comparison_case
    curves = TGConversionResult[r.cases[case] for r in runs]
    for setting in options.derivative_settings
        setting == options.primary_derivative && continue
        push!(
            profiles,
            fit_energy(:friedman, curves, options; case, derivative=setting, clocks),
        )
    end
    for delta in options.increments
        delta == options.protocol["primary_delta_alpha"] && continue
        push!(profiles, fit_energy(:advanced_vyazovkin, curves, options; case, delta))
    end
    refined = fit_energy(
        :advanced_vyazovkin,
        curves,
        options;
        case,
        subdivisions=options.protocol["refinement_subdivisions"],
    )
    for method in RAMP_METHODS
        push!(
            profiles,
            fit_energy(method, curves, options; case, ramp_rates=getfield.(runs, :rate)),
        )
    end
    foreach(require_estimated, profiles)
    require_estimated(refined)
    primary = only(
        filter(
            p ->
                p["case"] == case &&
                p["method"] == "advanced_vyazovkin" &&
                p["settings"]["delta_alpha"] == options.protocol["primary_delta_alpha"],
            profiles,
        ),
    )
    refinement_shift = maximum(abs, energies(refined) .- energies(primary))
    refinement_shift <= options.protocol["maximum_pp_refinement_shift_kJ_per_mol"] ||
        error("recorded-time energy refinement exceeds the declared tolerance")
    return Dict(
        "id" => cohort,
        "run_ids" => getfield.(runs, :id),
        "profiles" => profiles,
        "refined_primary" => refined,
        "max_refinement_shift_kJ_per_mol" => refinement_shift,
    )
end

function audit_public_tg_energy(output)
    protocol = TOML.parsefile(joinpath(TG_ROOT, ENERGY_PROTOCOL))
    conversion = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    derivative = TOML.parsefile(joinpath(TG_ROOT, FRIEDMAN_PROTOCOL))
    intervals = TOML.parsefile(joinpath(TG_ROOT, "config", "public_tg_vyazovkin.toml"))
    defaults = load_config(joinpath(TG_ROOT, "config", "analysis_defaults.toml"))
    options = (;
        protocol,
        advanced=defaults.isoconversional.advanced_vyazovkin,
        clock_bounds=derivative["clock_scale_alpha"],
        minimum_points=derivative["minimum_points_per_side"],
        primary_derivative=only(
            filter(s -> s["id"] == protocol["primary_derivative"], derivative["settings"])
        ),
        derivative_settings=derivative["settings"],
        increments=intervals["delta_alpha_grid"],
        omitted_run=derivative["omitted_polypropylene_run"],
    )
    controls = audit_controls(options)
    println(
        "Known-truth recovery gate passed. Estimating apparent energy for polypropylene."
    )
    runs = pp_runs(conversion, options)
    cohorts = [
        audit_pp_cohort("all_five", runs, conversion["comparison_case"], options),
        audit_pp_cohort(
            "omit_10",
            filter(r -> r.id != options.omitted_run, runs),
            conversion["comparison_case"],
            options,
        ),
    ]
    input_files = (
        ENERGY_PROTOCOL,
        TG_IMPORT_MANIFEST,
        TG_PROTOCOL,
        FRIEDMAN_PROTOCOL,
        "config/public_tg_vyazovkin.toml",
        "config/analysis_defaults.toml",
        "scripts/generate_tg_energy_controls.jl",
    )
    report = Dict(
        "schema_version" => 1,
        "protocol" => ENERGY_PROTOCOL,
        "alpha_grid" => protocol["alpha_grid"],
        "comparison_case" => conversion["comparison_case"],
        "known_truth_gate_passed" => true,
        "input_sha256" => Dict(
            path => bytes2hex(open(sha256, joinpath(TG_ROOT, path))) for
            path in input_files
        ),
        "controls" => controls,
        "polypropylene_cohorts" => cohorts,
        "polypropylene_runs" => [
            Dict(
                "id" => r.id,
                "path" => r.path,
                "sha256" => r.sha256,
                "full_record_OLS_K_per_min" => r.rate,
                "terminal_reference_mg" => r.terminal_mg,
            ) for r in runs
        ],
    )
    open(output, "w") do io
        return TOML.print(io, report; sorted=true)
    end
    println("Wrote $output")
    return report
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("Expected at most one report path")
    audit_public_tg_energy(
        isempty(ARGS) ? joinpath(TG_ROOT, ENERGY_REPORT) : abspath(only(ARGS))
    )
end
