tg_table_cell(value) = replace(string(value), '\n' => ' ', '\r' => ' ', '|' => "\\|")

function tg_table_row(io, values...)
    return println(io, "| ", join(tg_table_cell.(values), " | "), " |")
end

function tg_summary_value(mean, deviation)
    format = ICTACTGAudit.number
    isnothing(mean) && return "—"
    isnothing(deviation) && return format(mean)
    return "$(format(mean)) ± $(format(deviation))"
end

function render_tg_assessment(io, report, evidence_link)
    ICTACTGAudit.render_report(
        io,
        Dict("cohorts" => [report["assessment"]], "references" => report["references"]);
        generator="scripts/assess_tg.jl",
        evidence_link,
    )

    println(io, "## Thermal program roster\n")
    tg_table_row(
        io, "Program", "Role", "Kind", "Nominal rate (K/min)", "Nominal mass (mg)", "Runs"
    )
    tg_table_row(io, "---", "---", "---", "---:", "---:", "---")
    for program in report["programs"]
        tg_table_row(
            io,
            program["id"],
            replace(program["role"], '_' => ' '),
            program["kind"],
            get(program, "nominal_rate_K_per_min", "—"),
            get(program, "nominal_mass_mg", "—"),
            join(program["run_ids"], ", "),
        )
    end

    if !isempty(report["repeatability"])
        println(io, "\n## Repeatability and sample-mass checks\n")
        println(
            io,
            "Values with at least two observations include the sample standard deviation.\n",
        )
        tg_table_row(
            io,
            "Program",
            "α",
            "Supported crossings",
            "Crossing temperature (K)",
            "Retained mass (%)",
        )
        tg_table_row(io, "---", "---:", "---:", "---:", "---:")
        for group in report["repeatability"]
            residue = tg_summary_value(
                group["residue_mean_percent"], get(group, "residue_sd_percent", nothing)
            )
            for crossing in group["crossings"]
                temperature = tg_summary_value(
                    get(crossing, "temperature_mean_K", nothing),
                    get(crossing, "temperature_sd_K", nothing),
                )
                tg_table_row(
                    io,
                    group["program_id"],
                    ICTACTGAudit.number(crossing["alpha"]),
                    crossing["supported_repeats"],
                    temperature,
                    residue,
                )
            end
        end
    end

    if haskey(report, "prediction")
        prediction = report["prediction"]
        println(io, "\n## Frozen prediction on held-out programs\n")
        println(
            io,
            "One fit used $(length(prediction["training_ids"])) calibration runs. ",
            "That fit scored all $(length(prediction["held_out_ids"])) held-out runs.\n",
        )
        println(
            io,
            "Time-weighted RMSE is expressed as a percentage of the first mass used to ",
            "condition each prediction. The configuration declares the evaluation window.\n",
        )
        tg_table_row(
            io,
            "Program",
            "Run",
            "Evaluation RMSE (%)",
            "Full RMSE (%)",
            "Measured terminal mass (mg)",
            "Predicted terminal mass (mg)",
        )
        tg_table_row(io, "---", "---", "---:", "---:", "---:", "---:")
        format = ICTACTGAudit.number
        for score in prediction["hold_scores"]
            tg_table_row(
                io,
                score["program_id"],
                score["held_out_run"],
                format(score["evaluation"]["rmse_percent_first_mass"]),
                format(score["full"]["rmse_percent_first_mass"]),
                format(score["observed_final_mass_mg"]),
                format(score["predicted_final_mass_mg"]),
            )
        end
    end
    return nothing
end
