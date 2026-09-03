using CairoMakie
using IsoconversionalAnalysis
using Printf
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_REPORT = joinpath(PROJECT_ROOT, "docs", "preprocessing_audit.md")
const DEFAULT_MACHINE_REPORT = joinpath(
    PROJECT_ROOT, "docs", "audits", "m3_preprocessing_sensitivity.toml"
)
const DEFAULT_FIGURE_DIRECTORY = joinpath(PROJECT_ROOT, "docs", "figures", "m3")

function variant(
    original::PreprocessingConfig;
    profile=original.profile,
    invalid_row_policy=original.invalid_row_policy,
    segment_policy=original.segment_policy,
    rebase_time=original.rebase_time,
    reference_mass_method=original.reference_mass_method,
    reference_half_window_K=original.reference_half_window_K,
    smoothing_method=original.smoothing_method,
    smoothing_half_window_K=original.smoothing_half_window_K,
    local_polynomial_degree=original.local_polynomial_degree,
    derivative_method=original.derivative_method,
    heating_rate_method=original.heating_rate_method,
    monotonic_conversion_policy=original.monotonic_conversion_policy,
    analysis_conversion_range=original.analysis_conversion_range,
    reconstruction_rmse_tolerance=original.reconstruction_rmse_tolerance,
)
    return PreprocessingConfig(
        profile,
        invalid_row_policy,
        segment_policy,
        rebase_time,
        reference_mass_method,
        reference_half_window_K,
        smoothing_method,
        smoothing_half_window_K,
        local_polynomial_degree,
        derivative_method,
        heating_rate_method,
        monotonic_conversion_policy,
        analysis_conversion_range,
        reconstruction_rmse_tolerance,
    )
end

function with_preprocessing(config::AnalysisConfig, preprocessing::PreprocessingConfig)
    return AnalysisConfig(
        config.source_path,
        config.project,
        config.units,
        config.conversion,
        preprocessing,
        config.isoconversional,
        config.deconvolution,
        config.reaction_models,
    )
end

function ramp_maximum_celsius(experiment::Experiment, config::PreprocessingConfig)
    indices = findall(valid_row_mask(experiment))
    if config.segment_policy == :recorded_ramp && !isnothing(experiment.segment)
        target = experiment.kind == :dynamic ? 2 : 3
        recorded = filter(indices) do index
            segment = experiment.segment[index]
            return !ismissing(segment) && segment == target
        end
        !isempty(recorded) && (indices = recorded)
    end
    return maximum(experiment.temperature_K[indices]) - 273.15
end

function nearest_value(temperature_K, values, target_celsius)
    target_K = target_celsius + 273.15
    index = argmin(abs.(temperature_K .- target_K))
    return values[index]
end

function result_metrics(result::ProcessedExperiment)
    diagnostics = result.diagnostics
    row_count = length(result.alpha)
    roughness = sum(abs, diff(result.dalpha_dT_K_inv)) / max(row_count - 1, 1)
    return (
        reconstruction_rmse=diagnostics.reconstruction_rmse,
        reconstruction_max_error=diagnostics.reconstruction_max_abs_error,
        negative_derivative_fraction=diagnostics.negative_derivative_count / row_count,
        reversal_fraction=diagnostics.conversion_reversal_count / max(row_count - 1, 1),
        derivative_roughness=roughness,
        derivative_peak=maximum(result.dalpha_dT_K_inv),
    )
end

function median_metrics(results)
    metrics = result_metrics.(results)
    names = propertynames(first(metrics))
    values = map(names) do name
        return median(getproperty.(metrics, name))
    end
    return NamedTuple{names}(values)
end

function default_results(experiments, config)
    return Dict(
        experiment.id => preprocess(experiment, config) for experiment in experiments
    )
end

function candidate_audit(experiments, config, defaults)
    representative = filter(
        experiment -> experiment.composition.waste_tire == 0.5, experiments
    )
    base = config.preprocessing
    candidates = [
        variant(
            base;
            profile="raw_finite_difference",
            smoothing_method=:none,
            derivative_method=:finite_difference,
            heating_rate_method=:global_linear,
        ),
        variant(
            base;
            profile="smoothed_finite_difference_5K",
            derivative_method=:finite_difference,
        ),
        variant(base; profile="local_polynomial_2_5K", smoothing_half_window_K=2.5),
        base,
        variant(base; profile="local_polynomial_10K", smoothing_half_window_K=10.0),
    ]
    rows = NamedTuple[]
    result_sets = Dict{String,Vector{ProcessedExperiment}}()
    for candidate in candidates
        results = if candidate == base
            [defaults[experiment.id] for experiment in representative]
        else
            candidate_config = with_preprocessing(config, candidate)
            [preprocess(experiment, candidate_config) for experiment in representative]
        end
        result_sets[candidate.profile] = results
        push!(
            rows,
            (
                profile=candidate.profile,
                smoothing=String(candidate.smoothing_method),
                derivative=String(candidate.derivative_method),
                heating_rate=String(candidate.heating_rate_method),
                half_window_K=candidate.smoothing_half_window_K,
                metrics=median_metrics(results),
            ),
        )
    end
    return rows, result_sets
end

function reference_audit(experiments, config, defaults)
    linear = variant(
        config.preprocessing;
        profile="linear_reference_sensitivity",
        reference_mass_method=:linear_interpolation,
    )
    linear_config = with_preprocessing(config, linear)
    rows = NamedTuple[]
    for experiment in experiments
        default = defaults[experiment.id]
        comparison = preprocess(experiment, linear_config)
        push!(
            rows,
            (
                id=experiment.id,
                initial_mass_difference=comparison.diagnostics.initial_mass_percent -
                                        default.diagnostics.initial_mass_percent,
                final_mass_difference=comparison.diagnostics.final_mass_percent -
                                      default.diagnostics.final_mass_percent,
                alpha_350_difference=nearest_value(
                    comparison.temperature_K, comparison.alpha, 350.0
                ) - nearest_value(
                    default.temperature_K, default.alpha, 350.0
                ),
                alpha_500_difference=nearest_value(
                    comparison.temperature_K, comparison.alpha, 500.0
                ) - nearest_value(
                    default.temperature_K, default.alpha, 500.0
                ),
            ),
        )
    end
    return rows
end

function _temperature_at_alpha(result::ProcessedExperiment, target::Float64)
    index = argmin(abs.(result.analysis_alpha .- target))
    return result.temperature_K[index] - 273.15
end

function corrected_endpoint_audit(experiments, config, defaults, highest_common_celsius)
    profiles = [
        (name="initial_120_final_700", initial=120.0, final=700.0),
        (name="initial_150_final_650", initial=150.0, final=650.0),
        (name="primary_150_final_700", initial=150.0, final=700.0),
        (name="initial_150_final_750", initial=150.0, final=750.0),
        (name="initial_150_highest_common", initial=150.0, final=highest_common_celsius),
    ]
    checkpoints = [250.0, 350.0, 450.0, 550.0, 650.0]
    rows = NamedTuple[]
    for profile in profiles
        differences = Float64[]
        midpoint_shifts = Float64[]
        for experiment in experiments
            primary = defaults[experiment.id]
            result = if profile.name == "primary_150_final_700"
                primary
            else
                preprocess(
                    experiment,
                    config;
                    initial_temperature_celsius=profile.initial,
                    final_temperature_celsius=profile.final,
                )
            end
            for temperature in checkpoints
                profile.initial <= temperature <= min(profile.final, 700.0) || continue
                push!(
                    differences,
                    abs(
                        nearest_value(result.temperature_K, result.alpha, temperature) -
                        nearest_value(primary.temperature_K, primary.alpha, temperature),
                    ),
                )
            end
            push!(
                midpoint_shifts,
                _temperature_at_alpha(result, 0.5) - _temperature_at_alpha(primary, 0.5),
            )
        end
        push!(
            rows,
            (
                name=profile.name,
                initial_celsius=profile.initial,
                final_celsius=profile.final,
                median_absolute_alpha_difference=median(differences),
                maximum_absolute_alpha_difference=maximum(differences),
                median_temperature_at_alpha_half_shift_K=median(midpoint_shifts),
            ),
        )
    end
    return rows
end

function plot_qc(experiment, result, candidate_sets, output_path)
    set_theme!(Theme(; fontsize=15, linewidth=2))
    figure = Figure(; size=(1300, 900))
    finite = valid_row_mask(experiment)
    raw_temperature = experiment.temperature_K[finite] .- 273.15
    raw_mass = experiment.mass_percent[finite]
    retained_temperature = result.temperature_K .- 273.15

    mass_axis = Axis(
        figure[1, 1];
        xlabel="Temperature (°C)",
        ylabel="Mass (%)",
        title="Raw and retained mass",
    )
    lines!(mass_axis, raw_temperature, raw_mass; color=(:gray55, 0.65), label="finite raw")
    lines!(
        mass_axis,
        retained_temperature,
        result.mass_percent;
        color=:navy,
        label="150–700 °C",
    )
    vlines!(mass_axis, [150.0, 700.0]; color=:black, linestyle=:dash)
    scatter!(
        mass_axis,
        [150.0, 700.0],
        [result.diagnostics.initial_mass_percent, result.diagnostics.final_mass_percent];
        color=:darkorange,
        markersize=11,
        label="robust references",
    )
    axislegend(mass_axis; position=:rb)

    conversion_axis = Axis(
        figure[1, 2]; xlabel="Temperature (°C)", ylabel="Conversion, α", title="Conversion"
    )
    lines!(
        conversion_axis,
        retained_temperature,
        result.alpha;
        color=(:gray40, 0.6),
        label="raw",
    )
    lines!(
        conversion_axis,
        retained_temperature,
        result.analysis_alpha;
        color=:navy,
        label="local polynomial",
    )
    hlines!(conversion_axis, [0.0, 1.0]; color=:black, linestyle=:dot)
    axislegend(conversion_axis; position=:rb)

    rate_axis = Axis(
        figure[2, 1];
        xlabel="Temperature (°C)",
        ylabel="Heating rate (K/min)",
        title="Fitted local heating rate",
    )
    lines!(rate_axis, retained_temperature, result.heating_rate_K_per_min; color=:seagreen)
    hlines!(
        rate_axis,
        [experiment.nominal_heating_rate_K_per_min];
        color=:black,
        linestyle=:dash,
    )

    derivative_axis = Axis(
        figure[2, 2];
        xlabel="Temperature (°C)",
        ylabel="dα/dT (K⁻¹)",
        title="Derivative sensitivity",
    )
    palette = [:gray45, :darkorange, :dodgerblue3, :navy, :seagreen]
    for (color, profile) in zip(palette, sort(collect(keys(candidate_sets))))
        candidate = only(
            filter(item -> item.source_id == experiment.id, candidate_sets[profile])
        )
        lines!(
            derivative_axis,
            candidate.temperature_K .- 273.15,
            candidate.dalpha_dT_K_inv;
            color,
            label=profile,
        )
    end
    axislegend(derivative_axis; position=:rt, labelsize=10)

    reconstruction_axis = Axis(
        figure[3, 1];
        xlabel="Temperature (°C)",
        ylabel="Conversion, α",
        title="Integration reconstruction",
    )
    lines!(
        reconstruction_axis,
        retained_temperature,
        result.analysis_alpha;
        color=:navy,
        label="analysis α",
    )
    lines!(
        reconstruction_axis,
        retained_temperature,
        result.reconstructed_alpha;
        color=:darkorange,
        linestyle=:dash,
        label="∫(dα/dT)dT",
    )
    axislegend(reconstruction_axis; position=:rb)

    residual_axis = Axis(
        figure[3, 2];
        xlabel="Temperature (°C)",
        ylabel="Reconstruction residual",
        title="Reconstructed − analysis α",
    )
    lines!(
        residual_axis,
        retained_temperature,
        result.reconstruction_residual;
        color=:firebrick,
    )
    hlines!(residual_axis, [0.0]; color=:black, linestyle=:dot)

    Label(figure[0, :], "M3 preprocessing QC — $(experiment.id)"; fontsize=22, font=:bold)
    save(output_path, figure)
    return output_path
end

function plot_sensitivity(candidate_rows, endpoint_rows, reference_rows, output_path)
    figure = Figure(; size=(1300, 650))
    labels = replace.(getproperty.(candidate_rows, :profile), "_" => "\n")
    roughness = [row.metrics.derivative_roughness for row in candidate_rows]
    reconstruction = [row.metrics.reconstruction_rmse for row in candidate_rows]

    candidate_axis = Axis(
        figure[1, 1];
        xticks=(1:length(labels), labels),
        ylabel="Median metric",
        yscale=log10,
        title="Derivative roughness and reconstruction",
        xticklabelsize=10,
    )
    barplot!(candidate_axis, 1:length(labels), roughness; color=:navy, label="roughness")
    scatter!(
        candidate_axis,
        1:length(labels),
        reconstruction;
        color=:darkorange,
        markersize=13,
        label="RMSE",
    )
    axislegend(candidate_axis; position=:rt)

    endpoint_labels = replace.(getproperty.(endpoint_rows, :name), "_" => "\n")
    endpoint_axis = Axis(
        figure[1, 2];
        xticks=(1:length(endpoint_labels), endpoint_labels),
        ylabel="Absolute conversion difference",
        title="Endpoint sensitivity versus 150–700 °C",
        xticklabelsize=9,
    )
    barplot!(
        endpoint_axis,
        1:length(endpoint_labels),
        getproperty.(endpoint_rows, :median_absolute_alpha_difference);
        color=:seagreen,
        label="median",
    )
    scatter!(
        endpoint_axis,
        1:length(endpoint_labels),
        getproperty.(endpoint_rows, :maximum_absolute_alpha_difference);
        color=:firebrick,
        markersize=12,
        label="maximum",
    )
    axislegend(endpoint_axis; position=:rt)

    reference_axis = Axis(
        figure[1, 3];
        xlabel="Robust − interpolated reference mass (percentage points)",
        ylabel="Count",
        title="Reference estimator sensitivity",
    )
    differences = -getproperty.(reference_rows, :final_mass_difference)
    hist!(reference_axis, differences; bins=10, color=:dodgerblue3)

    Label(figure[0, :], "M3 sensitivity summary"; fontsize=22, font=:bold)
    save(output_path, figure)
    return output_path
end

function write_markdown(
    path,
    experiments,
    defaults,
    candidate_rows,
    endpoint_rows,
    reference_rows,
    highest_common_celsius,
)
    open(path, "w") do io
        println(io, "# M3 preprocessing audit")
        println(io)
        println(
            io,
            "Generated by `scripts/audit_preprocessing.jl` from the checksum-verified M2 catalog.",
        )
        println(io, "Raw `Experiment` arrays were not modified.")
        println(io)
        println(io, "## Accepted default")
        println(io)
        println(io, "The accepted `m3_default_v1` profile uses:")
        println(io)
        println(
            io,
            "- inclusive 150–700 °C selection after finite-row and recorded-ramp selection;",
        )
        println(
            io,
            "- robust local-linear reference masses evaluated exactly at 150 and 700 °C in ±2 K windows;",
        )
        println(io, "- conversion without clipping or monotonic repair;")
        println(
            io,
            "- a cubic local polynomial in time selected through a ±5 K physical window;",
        )
        println(
            io,
            "- local-polynomial `dα/dt` and local heating rate, with `dα/dT = (dα/dt)/(dT/dt)`;",
        )
        println(
            io,
            "- explicit warnings for invalid rows, missing segment metadata, conversion reversals, negative derivatives, and reconstruction failures.",
        )
        println(io)
        println(
            io,
            "The raw finite-difference profile remains available only as an audit comparison. No profile clips conversion or forces monotonicity.",
        )
        println(io)
        println(io, "## Acceptance checks")
        println(io)
        maximum_rmse = maximum(
            result.diagnostics.reconstruction_rmse for result in values(defaults)
        )
        println(
            io,
            "- Synthetic linear conversion: derivative and integrated reconstruction maximum errors are below `1e-8`.",
        )
        @printf(
            io,
            "- All %d dynamic runs preprocess successfully; maximum real-data reconstruction RMSE is `%.6g` (limit `0.01`).\n",
            length(experiments),
            maximum_rmse,
        )
        println(
            io,
            "- Source row indices and the complete configuration SHA-256 fingerprint are retained in every result.",
        )
        println(
            io,
            "- Ramp-and-hold experiments remain withheld and were not used to select this profile.",
        )
        println(io)
        println(io, "## Smoothing and derivative comparison")
        println(io)
        println(
            io,
            "Four 50 wt% tire runs spanning 5, 10, 15, and 20 K/min were used. Metrics below are medians across those runs.",
        )
        println(io)
        println(
            io,
            "| Profile | Window (K) | Smoothing | Derivative | Reconstruction RMSE | Negative derivative | Reversal fraction | Roughness | Peak (K⁻¹) |",
        )
        println(io, "|---|---:|---|---|---:|---:|---:|---:|---:|")
        for row in candidate_rows
            @printf(
                io,
                "| `%s` | %.1f | `%s` | `%s` | %.6g | %.4f | %.4f | %.6g | %.6g |\n",
                row.profile,
                row.half_window_K,
                row.smoothing,
                row.derivative,
                row.metrics.reconstruction_rmse,
                row.metrics.negative_derivative_fraction,
                row.metrics.reversal_fraction,
                row.metrics.derivative_roughness,
                row.metrics.derivative_peak,
            )
        end
        println(io)
        println(
            io,
            "The ±5 K local-polynomial profile is the balance point: it strongly reduces roughness without materially attenuating the median peak. The ±2.5 K and ±10 K profiles bracket the retained sensitivity range.",
        )
        println(io)
        println(io, "## Default results for all calibration runs")
        println(io)
        println(
            io,
            "| Experiment | Rows | Fitted β (K/min) | m150 (%) | m700 (%) | α range | Negative dα/dT | Reconstruction RMSE | Warnings |",
        )
        println(io, "|---|---:|---:|---:|---:|---|---:|---:|---:|")
        for experiment in experiments
            result = defaults[experiment.id]
            diagnostics = result.diagnostics
            @printf(
                io,
                "| `%s` | %d | %.4f | %.5f | %.5f | %.5f–%.5f | %.3f%% | %.6g | %d |\n",
                experiment.id,
                diagnostics.retained_row_count,
                diagnostics.measured_heating_rate_K_per_min,
                diagnostics.initial_mass_percent,
                diagnostics.final_mass_percent,
                diagnostics.alpha_minimum,
                diagnostics.alpha_maximum,
                100 * diagnostics.negative_derivative_count /
                    diagnostics.retained_row_count,
                diagnostics.reconstruction_rmse,
                length(diagnostics.warnings),
            )
        end
        println(io)
        println(
            io,
            "Warnings are diagnostics, not suppressed failures. The three 0 wt% runs lacking segment metadata explicitly record the all-finite fallback before temperature cropping.",
        )
        println(io)
        println(io, "## Reference-mass estimator sensitivity")
        println(io)
        @printf(
            io,
            "Across 20 runs, changing from robust local-linear references to shape-preserving linear interpolation changes the initial mass by median `%.6g` percentage points and the final mass by median `%.6g` percentage points in absolute value.\n",
            median(abs.(getproperty.(reference_rows, :initial_mass_difference))),
            median(abs.(getproperty.(reference_rows, :final_mass_difference))),
        )
        @printf(
            io,
            "The maximum absolute conversion change at 350 and 500 °C is `%.6g`. Robust local fitting is retained because it evaluates the requested temperature exactly while reducing single-sample sensitivity.\n",
            maximum(
                abs,
                vcat(
                    getproperty.(reference_rows, :alpha_350_difference),
                    getproperty.(reference_rows, :alpha_500_difference),
                ),
            ),
        )
        println(io)
        println(io, "## Conversion-endpoint sensitivity")
        println(io)
        @printf(
            io,
            "The highest common recorded ramp temperature is `%.5f °C`; no endpoint is extrapolated.\n",
            highest_common_celsius
        )
        println(io)
        println(
            io,
            "| Profile | Initial (°C) | Final (°C) | Median |Δα| | Maximum |Δα| | Median ΔT at α=0.5 (K) |",
        )
        println(io, "|---|---:|---:|---:|---:|---:|")
        for row in endpoint_rows
            @printf(
                io,
                "| `%s` | %.3f | %.5f | %.6g | %.6g | %.5f |\n",
                row.name,
                row.initial_celsius,
                row.final_celsius,
                row.median_absolute_alpha_difference,
                row.maximum_absolute_alpha_difference,
                row.median_temperature_at_alpha_half_shift_K,
            )
        end
        println(io)
        println(
            io,
            "The 150–700 °C coordinate remains primary. M4 must propagate the alternative endpoint profiles as a preprocessing sensitivity rather than treating their spread as replicate uncertainty.",
        )
        println(io)
        println(io, "## Figures")
        println(io)
        println(
            io,
            "- [`figures/m3/preprocessing_qc.svg`](figures/m3/preprocessing_qc.svg): retained interval, reference masses, conversion, heating rate, derivative candidates, and integration residual.",
        )
        println(
            io,
            "- [`figures/m3/preprocessing_sensitivity.svg`](figures/m3/preprocessing_sensitivity.svg): candidate, endpoint, and reference-estimator sensitivity summaries.",
        )
        println(io)
        println(io, "## Legacy deviations")
        println(io)
        return println(
            io,
            "The rewrite intentionally does not reproduce `cleandata.m`'s first/last retained-sample normalization, strict endpoint exclusion, derivative stencil on irregular coordinates, or unexplained low-pass cutoff. These are documented comparison artifacts, not compatibility requirements.",
        )
    end
    return path
end

function toml_string(value::AbstractString)
    return "\"" * replace(value, "\\" => "\\\\", "\"" => "\\\"") * "\""
end

function write_machine_report(
    path,
    experiments,
    defaults,
    candidate_rows,
    endpoint_rows,
    reference_rows,
    highest_common_celsius,
)
    open(path, "w") do io
        println(io, "schema_version = 1")
        println(io, "generated_by = \"scripts/audit_preprocessing.jl\"")
        println(io, "default_profile = \"m3_default_v1\"")
        @printf(io, "highest_common_temperature_celsius = %.10f\n", highest_common_celsius)
        println(io)
        for row in candidate_rows
            println(io, "[[candidate]]")
            println(io, "profile = ", toml_string(row.profile))
            println(io, "smoothing = ", toml_string(row.smoothing))
            println(io, "derivative = ", toml_string(row.derivative))
            println(io, "heating_rate = ", toml_string(row.heating_rate))
            @printf(io, "half_window_kelvin = %.10f\n", row.half_window_K)
            for name in propertynames(row.metrics)
                @printf(io, "%s = %.16g\n", String(name), getproperty(row.metrics, name))
            end
            println(io)
        end
        for row in endpoint_rows
            println(io, "[[endpoint]]")
            println(io, "profile = ", toml_string(row.name))
            @printf(io, "initial_temperature_celsius = %.10f\n", row.initial_celsius)
            @printf(io, "final_temperature_celsius = %.10f\n", row.final_celsius)
            @printf(
                io,
                "median_absolute_alpha_difference = %.16g\n",
                row.median_absolute_alpha_difference
            )
            @printf(
                io,
                "maximum_absolute_alpha_difference = %.16g\n",
                row.maximum_absolute_alpha_difference
            )
            @printf(
                io,
                "median_temperature_at_alpha_half_shift_kelvin = %.16g\n",
                row.median_temperature_at_alpha_half_shift_K
            )
            println(io)
        end
        for row in reference_rows
            println(io, "[[reference_estimator]]")
            println(io, "experiment_id = ", toml_string(row.id))
            @printf(
                io, "initial_mass_difference_percent = %.16g\n", row.initial_mass_difference
            )
            @printf(
                io, "final_mass_difference_percent = %.16g\n", row.final_mass_difference
            )
            @printf(io, "alpha_350_difference = %.16g\n", row.alpha_350_difference)
            @printf(io, "alpha_500_difference = %.16g\n", row.alpha_500_difference)
            println(io)
        end
        for experiment in experiments
            result = defaults[experiment.id]
            diagnostics = result.diagnostics
            println(io, "[[experiment]]")
            println(io, "id = ", toml_string(experiment.id))
            println(io, "config_fingerprint = ", toml_string(result.config_fingerprint))
            println(io, "retained_rows = ", diagnostics.retained_row_count)
            @printf(io, "initial_mass_percent = %.16g\n", diagnostics.initial_mass_percent)
            @printf(io, "final_mass_percent = %.16g\n", diagnostics.final_mass_percent)
            @printf(
                io,
                "fitted_heating_rate_kelvin_per_minute = %.16g\n",
                diagnostics.measured_heating_rate_K_per_min
            )
            @printf(io, "alpha_minimum = %.16g\n", diagnostics.alpha_minimum)
            @printf(io, "alpha_maximum = %.16g\n", diagnostics.alpha_maximum)
            println(
                io,
                "outside_unit_interval_count = ",
                diagnostics.alpha_outside_unit_interval_count,
            )
            println(
                io, "conversion_reversal_count = ", diagnostics.conversion_reversal_count
            )
            println(
                io, "negative_derivative_count = ", diagnostics.negative_derivative_count
            )
            @printf(io, "reconstruction_rmse = %.16g\n", diagnostics.reconstruction_rmse)
            @printf(
                io,
                "reconstruction_max_abs_error = %.16g\n",
                diagnostics.reconstruction_max_abs_error
            )
            println(io, "warnings = [", join(toml_string.(diagnostics.warnings), ", "), "]")
            println(io)
        end
    end
    return path
end

function main()
    config = load_config(joinpath(PROJECT_ROOT, "config", "analysis_defaults.toml"))
    catalog = load_dataset_catalog(joinpath(PROJECT_ROOT, "config", "datasets.toml"))
    experiments = load_experiments(catalog; roles=[:calibration])
    defaults = default_results(experiments, config)
    highest_common_celsius = minimum(
        ramp_maximum_celsius(experiment, config.preprocessing) for experiment in experiments
    )

    candidate_rows, candidate_sets = candidate_audit(experiments, config, defaults)
    reference_rows = reference_audit(experiments, config, defaults)
    endpoint_rows = corrected_endpoint_audit(
        experiments, config, defaults, highest_common_celsius
    )

    mkpath(dirname(DEFAULT_REPORT))
    mkpath(dirname(DEFAULT_MACHINE_REPORT))
    mkpath(DEFAULT_FIGURE_DIRECTORY)
    write_markdown(
        DEFAULT_REPORT,
        experiments,
        defaults,
        candidate_rows,
        endpoint_rows,
        reference_rows,
        highest_common_celsius,
    )
    write_machine_report(
        DEFAULT_MACHINE_REPORT,
        experiments,
        defaults,
        candidate_rows,
        endpoint_rows,
        reference_rows,
        highest_common_celsius,
    )

    representative_experiment = only(
        filter(experiments) do experiment
            return experiment.id == "dynamic_wt050_rate10"
        end,
    )
    plot_qc(
        representative_experiment,
        defaults[representative_experiment.id],
        candidate_sets,
        joinpath(DEFAULT_FIGURE_DIRECTORY, "preprocessing_qc.svg"),
    )
    plot_sensitivity(
        candidate_rows,
        endpoint_rows,
        reference_rows,
        joinpath(DEFAULT_FIGURE_DIRECTORY, "preprocessing_sensitivity.svg"),
    )
    println("Wrote $(relpath(DEFAULT_REPORT, PROJECT_ROOT))")
    println("Wrote $(relpath(DEFAULT_MACHINE_REPORT, PROJECT_ROOT))")
    return println("Wrote $(relpath(DEFAULT_FIGURE_DIRECTORY, PROJECT_ROOT))/")
end

main()
