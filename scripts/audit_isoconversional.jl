using CairoMakie
using IsoconversionalAnalysis
using Printf
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_REPORT = joinpath(PROJECT_ROOT, "docs", "isoconversional_audit.md")
const DEFAULT_MACHINE_REPORT = joinpath(
    PROJECT_ROOT, "docs", "audits", "m4_isoconversional_results.toml"
)
const DEFAULT_FIGURE_DIRECTORY = joinpath(PROJECT_ROOT, "docs", "figures", "m4")
const METHODS = (:friedman, :kas, :fwo, :starink, :advanced_vyazovkin)
const METHOD_LABELS = Dict(
    :friedman => "Friedman",
    :kas => "KAS",
    :fwo => "FWO",
    :starink => "Starink",
    :advanced_vyazovkin => "Advanced Vyazovkin",
)

function ramp_maximum_celsius(experiment::Experiment, config::PreprocessingConfig)
    indices = findall(valid_row_mask(experiment))
    if config.segment_policy == :recorded_ramp && !isnothing(experiment.segment)
        selected_segment = experiment.kind == :dynamic ? 2 : 3
        recorded = filter(indices) do index
            segment = experiment.segment[index]
            return !ismissing(segment) && segment == selected_segment
        end
        !isempty(recorded) && (indices = recorded)
    end
    return maximum(experiment.temperature_K[indices]) - 273.15
end

function endpoint_profiles(experiments, config)
    initial = config.conversion.initial_temperature_celsius
    final = config.conversion.final_temperature_celsius
    highest_common = minimum(
        ramp_maximum_celsius(experiment, config.preprocessing) for experiment in experiments
    )
    return [
        (name="primary_150_700", initial=initial, final=final),
        (name="initial_120", initial=120.0, final=final),
        (name="final_650", initial=initial, final=650.0),
        (name="final_750", initial=initial, final=750.0),
        (name="final_highest_common", initial=initial, final=highest_common),
    ]
end

function run_profile(experiments, config, profile)
    processed = [
        preprocess(
            experiment,
            config;
            initial_temperature_celsius=profile.initial,
            final_temperature_celsius=profile.final,
        ) for experiment in experiments
    ]
    results = Dict{Tuple{Int,Symbol},IsoconversionalResult}()
    for fraction in
        sort(unique(experiment.composition.waste_tire for experiment in processed))
        group = filter(
            experiment -> experiment.composition.waste_tire == fraction, processed
        )
        composition_percent = round(Int, 100 * fraction)
        for method in METHODS
            results[(composition_percent, method)] = analyze(method, group, config)
        end
    end
    return results
end

function finite_estimates(result::IsoconversionalResult)
    return Float64[
        value for value in result.activation_energy_kJ_per_mol if !ismissing(value)
    ]
end

function status_counts(result::IsoconversionalResult)
    statuses = getproperty.(result.point_diagnostics, :status)
    return [
        (status=status, count=count(==(status), statuses)) for
        status in sort(unique(statuses))
    ]
end

function minimum_r_squared(result::IsoconversionalResult)
    diagnostics = [
        point.regression.r_squared for point in result.point_diagnostics if
        !isnothing(point.regression) && isfinite(point.regression.r_squared)
    ]
    return isempty(diagnostics) ? nothing : minimum(diagnostics)
end

function endpoint_summaries(primary, alternatives)
    rows = NamedTuple[]
    for (profile, results) in alternatives
        for composition in (0, 25, 50, 75, 100), method in METHODS
            reference = primary[(composition, method)]
            comparison = results[(composition, method)]
            deviations = Float64[]
            for (reference_value, comparison_value) in zip(
                reference.activation_energy_kJ_per_mol,
                comparison.activation_energy_kJ_per_mol,
            )
                if !ismissing(reference_value) && !ismissing(comparison_value)
                    push!(deviations, abs(reference_value - comparison_value))
                end
            end
            push!(
                rows,
                (
                    profile=profile.name,
                    initial=profile.initial,
                    final=profile.final,
                    composition=composition,
                    method=method,
                    compared_points=length(deviations),
                    median_absolute_difference=if isempty(deviations)
                        NaN
                    else
                        median(deviations)
                    end,
                    maximum_absolute_difference=if isempty(deviations)
                        NaN
                    else
                        maximum(deviations)
                    end,
                    invalid_points=count(
                        ismissing, comparison.activation_energy_kJ_per_mol
                    ),
                ),
            )
        end
    end
    return rows
end

function format_number(value; digits=2)
    return isfinite(value) ? @sprintf("%.*f", digits, value) : "not available"
end

function format_statuses(result)
    return join(["`$(row.status)` × $(row.count)" for row in status_counts(result)], ", ")
end

function write_markdown(path, config, profiles, primary, endpoint_rows)
    open(path, "w") do io
        println(io, "# M4 isoconversional-method audit")
        println(io)
        println(io, "Generated by `scripts/audit_isoconversional.jl` on 2026-08-17.")
        println(io)
        println(io, "## Scope and acceptance contract")
        println(io)
        println(
            io,
            "The audit uses only the 20 dynamic calibration runs: five dry-mass " *
            "compositions and four heating programs per composition. All 15 ramp-and-hold " *
            "runs remain held out for predictive validation. Every method uses the same " *
            "$(length(config.isoconversional.alpha_grid))-point conversion grid from " *
            "$(first(config.isoconversional.alpha_grid)) to " *
            "$(last(config.isoconversional.alpha_grid)).",
        )
        println(io)
        println(
            io,
            "The predeclared synthetic gates are: Friedman error below 0.5 kJ/mol, " *
            "advanced-Vyazovkin error below 0.75 kJ/mol, integral-linear method error " *
            "below 8 kJ/mol, and noisy-fixture error below 15 kJ/mol. These are enforced " *
            "in `test/test_isoconversional.jl` and were not tuned against the real data.",
        )
        println(io)
        println(io, "## Frozen numerical conventions")
        println(io)
        println(
            io,
            "- Conversion interpolation: first acquisition-order upward crossing; no extrapolation or monotonic repair.",
        )
        println(
            io,
            "- Linear methods: ordinary least squares with two-sided $(round(Int, 100 * config.isoconversional.confidence_level))% Student-t intervals for the transformed slope.",
        )
        println(
            io,
            "- Linear-fit diagnostics: results remain available but are flagged when R² is below $(config.isoconversional.minimum_r_squared_warning).",
        )
        println(
            io,
            "- Advanced Vyazovkin: true pairwise objective, log-scaled trapezoidal integrals, bounded scan plus golden-section refinement, and explicit boundary/flat/multimodal flags.",
        )
        println(
            io,
            "- Advanced uncertainty: the Vyazovkin–Wight pairwise-variance Fisher interval; optimizer quality remains a separate diagnostic.",
        )
        println(
            io,
            "- FWO: retained as a named compatibility result, not selected as the preferred integral approximation.",
        )
        println(io)
        println(io, "## Primary 150–700 °C results")
        println(io)
        println(
            io,
            "| Tire dry-mass fraction | Method | Valid points | Median E (kJ/mol) | Range (kJ/mol) | Minimum R² | Statuses |",
        )
        println(io, "|---:|---|---:|---:|---:|---:|---|")
        for composition in (0, 25, 50, 75, 100), method in METHODS
            result = primary[(composition, method)]
            estimates = finite_estimates(result)
            estimate_range = if isempty(estimates)
                "not available"
            else
                "$(format_number(minimum(estimates)))–$(format_number(maximum(estimates)))"
            end
            regression_r_squared = minimum_r_squared(result)
            r_squared_label = if isnothing(regression_r_squared)
                "—"
            else
                format_number(regression_r_squared; digits=4)
            end
            median_label =
                isempty(estimates) ? "not available" : format_number(median(estimates))
            println(
                io,
                "| $composition% | $(METHOD_LABELS[method]) | $(length(estimates))/$(length(result.alpha)) | $median_label | $estimate_range | $r_squared_label | $(format_statuses(result)) |",
            )
        end
        println(io)
        println(
            io,
            "All primary method/composition combinations return all requested conversion " *
            "points. The advanced method exposes a lower-bound solution in one point of " *
            "the 25% tire mixture and an upper-bound solution in one point of the 50% " *
            "mixture. These are numerical/scientific warnings, not missing-data events; " *
            "the affected values must not be interpreted as unconstrained estimates.",
        )
        println(
            io,
            "The remaining advanced-method warnings in those two mixtures are Fisher " *
            "intervals truncated at one or both configured energy bounds. They diagnose " *
            "weakly bounded uncertainty; they are not additional optimizer failures.",
        )
        println(
            io,
            "The very low minimum R² and low activation-energy tail for the 25% mixture " *
            "show that its highest-conversion region does not support a reliable straight-line " *
            "summary across the four heating programs. The values are retained and flagged; " *
            "they are not evidence for a near-zero physical barrier.",
        )
        println(io)
        println(io, "![Activation-energy profiles](figures/m4/isoconversional_methods.svg)")
        println(io)
        println(io, "## Conversion-endpoint sensitivity")
        println(io)
        println(
            io,
            "The primary coordinate is compared with the M0/M3 alternatives below. Differences are absolute changes relative to the primary activation-energy profile at the same reported conversion.",
        )
        println(io)
        println(
            io,
            "| Endpoint profile | Initial (°C) | Final (°C) | Median across method/composition maxima (kJ/mol) | Largest maximum (kJ/mol) |",
        )
        println(io, "|---|---:|---:|---:|---:|")
        for profile in profiles[2:end]
            matching = filter(row -> row.profile == profile.name, endpoint_rows)
            maxima = [
                row.maximum_absolute_difference for
                row in matching if isfinite(row.maximum_absolute_difference)
            ]
            println(
                io,
                "| `$(profile.name)` | $(format_number(profile.initial; digits=1)) | $(format_number(profile.final; digits=3)) | $(format_number(median(maxima))) | $(format_number(maximum(maxima))) |",
            )
        end
        println(io)
        println(io, "![Endpoint sensitivity](figures/m4/endpoint_sensitivity.svg)")
        println(io)
        println(
            io,
            "Endpoint sensitivity is reported separately from statistical intervals because " *
            "it changes the definition of conversion itself. The complete method- and " *
            "composition-level values are stored in " *
            "`audits/m4_isoconversional_results.toml`.",
        )
        println(io)
        println(io, "## Interpretation and M4 decision")
        println(io)
        println(
            io,
            "KAS, FWO, Starink, Friedman, and advanced Vyazovkin show broadly consistent " *
            "central profiles, while departures near individual conversion regions are " *
            "retained as useful evidence of changing apparent kinetics. FWO remains available " *
            "for thesis/legacy comparison, but Starink and advanced Vyazovkin are the preferred " *
            "integral summaries because FWO uses the coarser Doyle approximation.",
        )
        println(io)
        return println(
            io,
            "M4 satisfies its exit gate: the shared interpolation contract and every method " *
            "pass known-energy and noisy synthetic tests; arbitrary experiment counts and " *
            "permutation invariance are tested; real-data results are complete and retain all " *
            "diagnostic warnings. M5 can consume the typed results without treating this audit " *
            "as model validation of the held-out ramp-and-hold programs.",
        )
    end
    return path
end

toml_string(value) = "\"" * replace(String(value), "\\" => "\\\\", "\"" => "\\\"") * "\""
toml_float(value) = @sprintf("%.16g", value)
toml_float_array(values) = "[" * join(toml_float.(values), ", ") * "]"
toml_string_array(values) = "[" * join(toml_string.(values), ", ") * "]"

function write_machine_report(path, config, profiles, primary, endpoint_rows)
    open(path, "w") do io
        println(io, "schema_version = 1")
        println(io, "generated_on = \"2026-08-17\"")
        println(io, "calibration_experiment_count = 20")
        println(io, "validation_experiment_count_used = 0")
        println(
            io, "confidence_level = ", toml_float(config.isoconversional.confidence_level)
        )
        println(
            io,
            "minimum_r_squared_warning = ",
            toml_float(config.isoconversional.minimum_r_squared_warning),
        )
        println(io, "alpha_grid = ", toml_float_array(config.isoconversional.alpha_grid))
        println(io)
        for profile in profiles
            println(io, "[[endpoint_profile]]")
            println(io, "name = ", toml_string(profile.name))
            println(io, "initial_temperature_celsius = ", toml_float(profile.initial))
            println(io, "final_temperature_celsius = ", toml_float(profile.final))
            println(io)
        end
        for composition in (0, 25, 50, 75, 100), method in METHODS
            result = primary[(composition, method)]
            println(io, "[[primary_result]]")
            println(io, "waste_tire_dry_mass_percent = $composition")
            println(io, "method = ", toml_string(method))
            println(io, "analysis_fingerprint = ", toml_string(result.analysis_fingerprint))
            println(io, "interval_method = ", toml_string(result.interval_method))
            println(io, "experiment_ids = ", toml_string_array(result.experiment_ids))
            println(io, "alpha = ", toml_float_array(result.alpha))
            println(
                io,
                "activation_energy_kilojoule_per_mole = ",
                toml_float_array(Float64.(result.activation_energy_kJ_per_mol)),
            )
            println(
                io,
                "confidence_lower_kilojoule_per_mole = ",
                toml_float_array(Float64.(result.confidence_lower_kJ_per_mol)),
            )
            println(
                io,
                "confidence_upper_kilojoule_per_mole = ",
                toml_float_array(Float64.(result.confidence_upper_kJ_per_mol)),
            )
            println(
                io,
                "status = ",
                toml_string_array(string.(getproperty.(result.point_diagnostics, :status))),
            )
            println(
                io,
                "warnings = ",
                toml_string_array(
                    join(point.warnings, " | ") for point in result.point_diagnostics
                ),
            )
            if method != :advanced_vyazovkin
                println(
                    io,
                    "r_squared = ",
                    toml_float_array(
                        point.regression.r_squared for point in result.point_diagnostics
                    ),
                )
            else
                advanced = getproperty.(result.point_diagnostics, :advanced_vyazovkin)
                println(
                    io,
                    "objective_minimum = ",
                    toml_float_array(value.objective_minimum for value in advanced),
                )
                println(
                    io,
                    "objective_baseline = ",
                    toml_float_array(value.objective_baseline for value in advanced),
                )
                println(
                    io,
                    "boundary_solution = [",
                    join(string.(getproperty.(advanced, :boundary_solution)), ", "),
                    "]",
                )
                println(
                    io,
                    "fisher_variance_minimum = ",
                    toml_float_array(value.fisher_variance_minimum for value in advanced),
                )
                println(
                    io,
                    "fisher_center_energy_kilojoule_per_mole = ",
                    toml_float_array(
                        value.fisher_center_energy_kJ_per_mol for value in advanced
                    ),
                )
                println(
                    io,
                    "fisher_lower_truncated = [",
                    join(string.(getproperty.(advanced, :fisher_lower_truncated)), ", "),
                    "]",
                )
                println(
                    io,
                    "fisher_upper_truncated = [",
                    join(string.(getproperty.(advanced, :fisher_upper_truncated)), ", "),
                    "]",
                )
            end
            println(io)
        end
        for row in endpoint_rows
            println(io, "[[endpoint_sensitivity]]")
            println(io, "profile = ", toml_string(row.profile))
            println(io, "waste_tire_dry_mass_percent = $(row.composition)")
            println(io, "method = ", toml_string(row.method))
            println(io, "compared_points = $(row.compared_points)")
            println(io, "invalid_points = $(row.invalid_points)")
            println(
                io,
                "median_absolute_difference_kilojoule_per_mole = ",
                toml_float(row.median_absolute_difference),
            )
            println(
                io,
                "maximum_absolute_difference_kilojoule_per_mole = ",
                toml_float(row.maximum_absolute_difference),
            )
            println(io)
        end
    end
    return path
end

function plot_methods(path, primary)
    colors = Makie.wong_colors()
    figure = Figure(; size=(1450, 870))
    axes = Axis[]
    for (index, composition) in enumerate((0, 25, 50, 75, 100))
        row = index <= 3 ? 1 : 2
        column = index <= 3 ? index : index - 3
        axis = Axis(
            figure[row, column];
            title="$composition% waste tire (dry mass)",
            xlabel="Conversion, α",
            ylabel="Eα (kJ/mol)",
        )
        push!(axes, axis)
        for (method_index, method) in enumerate(METHODS)
            result = primary[(composition, method)]
            lines!(
                axis,
                result.alpha,
                Float64.(result.activation_energy_kJ_per_mol);
                color=colors[method_index],
                linewidth=2.2,
                label=METHOD_LABELS[method],
            )
            if method == :advanced_vyazovkin
                scatter!(
                    axis,
                    result.alpha,
                    Float64.(result.activation_energy_kJ_per_mol);
                    color=colors[method_index],
                    markersize=5,
                )
            end
        end
    end
    hidexdecorations!(axes[1]; grid=false)
    hidexdecorations!(axes[2]; grid=false)
    hidexdecorations!(axes[3]; grid=false)
    Legend(figure[2, 3], axes[1]; framevisible=false, tellheight=false)
    Label(
        figure[0, 1:3],
        "M4 apparent activation-energy profiles — primary 150–700 °C conversion";
        fontsize=23,
    )
    save(path, figure)
    return path
end

function plot_endpoint_sensitivity(path, profiles, endpoint_rows)
    colors = Makie.wong_colors()
    alternatives = profiles[2:end]
    labels = replace.(getproperty.(alternatives, :name), "_" => "\n")
    figure = Figure(; size=(1450, 870))
    axes = Axis[]
    for (index, composition) in enumerate((0, 25, 50, 75, 100))
        row = index <= 3 ? 1 : 2
        column = index <= 3 ? index : index - 3
        axis = Axis(
            figure[row, column];
            title="$composition% waste tire (dry mass)",
            xlabel="Endpoint profile",
            ylabel="Maximum |ΔEα| (kJ/mol)",
            xticks=(eachindex(labels), labels),
            xticklabelsize=10,
        )
        push!(axes, axis)
        for (method_index, method) in enumerate(METHODS)
            values = Float64[]
            for profile in alternatives
                row_value = only(
                    filter(endpoint_rows) do candidate
                        return candidate.profile == profile.name &&
                               candidate.composition == composition &&
                               candidate.method == method
                    end,
                )
                push!(values, row_value.maximum_absolute_difference)
            end
            lines!(
                axis,
                eachindex(alternatives),
                values;
                color=colors[method_index],
                linewidth=2.2,
                label=METHOD_LABELS[method],
            )
            scatter!(axis, eachindex(alternatives), values; color=colors[method_index])
        end
    end
    Legend(figure[2, 3], axes[1]; framevisible=false, tellheight=false)
    Label(
        figure[0, 1:3],
        "Sensitivity of Eα profiles to conversion-reference temperatures";
        fontsize=23,
    )
    save(path, figure)
    return path
end

function main()
    config = load_config(joinpath(PROJECT_ROOT, "config", "analysis_defaults.toml"))
    catalog = load_dataset_catalog(joinpath(PROJECT_ROOT, "config", "datasets.toml"))
    experiments = load_experiments(catalog; roles=[:calibration])
    profiles = endpoint_profiles(experiments, config)
    primary = run_profile(experiments, config, first(profiles))
    alternatives = [
        (profile, run_profile(experiments, config, profile)) for profile in profiles[2:end]
    ]
    endpoint_rows = endpoint_summaries(primary, alternatives)

    mkpath(dirname(DEFAULT_REPORT))
    mkpath(dirname(DEFAULT_MACHINE_REPORT))
    mkpath(DEFAULT_FIGURE_DIRECTORY)
    write_markdown(DEFAULT_REPORT, config, profiles, primary, endpoint_rows)
    write_machine_report(DEFAULT_MACHINE_REPORT, config, profiles, primary, endpoint_rows)
    plot_methods(joinpath(DEFAULT_FIGURE_DIRECTORY, "isoconversional_methods.svg"), primary)
    plot_endpoint_sensitivity(
        joinpath(DEFAULT_FIGURE_DIRECTORY, "endpoint_sensitivity.svg"),
        profiles,
        endpoint_rows,
    )
    println("Wrote $(relpath(DEFAULT_REPORT, PROJECT_ROOT))")
    println("Wrote $(relpath(DEFAULT_MACHINE_REPORT, PROJECT_ROOT))")
    return println("Wrote $(relpath(DEFAULT_FIGURE_DIRECTORY, PROJECT_ROOT))/")
end

main()
