using CairoMakie
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const MATERIALS = (
    ("Constant_E", "Constant energy"),
    ("Two_Steps", "Two-step"),
    ("Variable_E", "Variable energy"),
    ("polypropylene", "Polypropylene"),
)
const COLORS = ("#27638A", "#B64B25", "#367D62", "#5F5089")
const MARKERS = (:circle, :rect, :diamond, :utriangle)

central(points) = filter(point -> 0.10 <= point["alpha"] <= 0.90, points)

function comparison_case(report, material; field="comparison_case")
    group = only(filter(group -> group["material_id"] == material, report["groups"]))
    return only(filter(case -> case["id"] == report[field], group["cases"]))
end

function method_figure(friedman, vyazovkin, ramp)
    fig = Figure(; size=(1240, 835), fontsize=15, backgroundcolor=:white)
    rate = Axis(
        fig[1, 1];
        title="Friedman derivative sensitivity",
        xlabel="Conversion, α",
        ylabel="Largest rate spread (%)",
        xticks=0.1:0.2:0.9,
    )
    for (i, (material, label)) in enumerate(MATERIALS)
        points = central(comparison_case(friedman, material)["points"])
        values = [100 * expm1(get(p, "max_within_run_log_rate_span", NaN)) for p in points]
        scatterlines!(
            rate,
            [p["alpha"] for p in points],
            values;
            color=COLORS[i],
            marker=MARKERS[i],
            markersize=8,
            label,
        )
    end
    Legend(fig[2, 1], rate; framevisible=false, labelsize=13, nbanks=2, tellwidth=false)

    alpha = vyazovkin["alpha_grid"]
    widths = vyazovkin["delta_alpha_grid"]
    cases = [comparison_case(vyazovkin, material) for (material, _) in MATERIALS]
    counts = [
        count(w -> a in w["common_supported_alpha"], c["windows"]) for
        a in alpha, c in cases
    ]
    ticks = [1, 4, 6, 8, 10, 12, length(alpha)]
    intervals = Axis(
        fig[1, 2];
        title="Advanced Vyazovkin interval support",
        xlabel="Tested conversion, α",
        xticks=(ticks, string.(alpha[ticks])),
        yticks=(1:length(MATERIALS), [label for (_, label) in MATERIALS]),
        yreversed=true,
        xgridvisible=false,
        ygridvisible=false,
    )
    palette = ["#E9ECEF", "#C7DDE1", "#81AFBB", "#2B6981"]
    heatmap!(
        intervals,
        1:length(alpha),
        1:length(MATERIALS),
        counts;
        colormap=cgrad(palette; categorical=true),
        colorrange=(0, length(widths)),
    )
    for i in axes(counts, 1), j in axes(counts, 2)
        n = counts[i, j]
        text!(
            intervals,
            i,
            j;
            text=string(n),
            align=(:center, :center),
            fontsize=11,
            color=n == length(widths) ? :white : :black,
        )
    end
    Label(
        fig[2, 2],
        "Number of supported Δα widths (0.01, 0.02, 0.05)";
        fontsize=13,
        tellwidth=false,
    )

    pp = comparison_case(ramp, "polypropylene"; field="baseline_case")
    start = only(filter(start -> start["id"] == "selected_first_row", pp["program_starts"]))
    points = filter(p -> p["all_prefixes_evaluated"], central(start["points"]))
    prefix_alpha = [p["alpha"] for p in points]
    path_axis = Axis(
        fig[3, 1];
        title="PP recorded versus affine ramp",
        xlabel="Conversion, α",
        ylabel="Largest |log(integral ratio)|",
        xticks=0.1:0.2:0.9,
    )
    for (i, energy) in enumerate(ramp["probe_energies_kJ_per_mol"])
        label = "$(Int(energy)) kJ/mol probe"
        values = [
            p["probes"][i]["max_absolute_measured_vs_affine_log_ratio"] for p in points
        ]
        scatterlines!(
            path_axis,
            prefix_alpha,
            values;
            color=COLORS[i],
            marker=MARKERS[i],
            markersize=8,
            label,
        )
    end
    Legend(
        fig[4, 1], path_axis; framevisible=false, labelsize=13, nbanks=2, tellwidth=false
    )

    shape_axis = Axis(
        fig[3, 2];
        title="PP formula shape at 120 kJ/mol",
        xlabel="Conversion, α",
        ylabel="Across-run log residual span",
        xticks=0.1:0.2:0.9,
    )
    energy_index = only(findall(==(120), ramp["probe_energies_kJ_per_mol"]))
    for (i, label) in enumerate(("KAS", "FWO", "Starink"))
        values = [
            p["probes"][energy_index]["approximation_log_residual_span"][lowercase(label)]
            for p in points
        ]
        scatterlines!(
            shape_axis,
            prefix_alpha,
            values;
            color=COLORS[i],
            marker=MARKERS[i],
            markersize=8,
            label,
        )
    end
    Legend(
        fig[4, 2],
        shape_axis;
        framevisible=false,
        labelsize=13,
        orientation=:horizontal,
        tellwidth=false,
    )
    rowsize!(fig.layout, 2, Fixed(45))
    rowsize!(fig.layout, 4, Fixed(45))
    rowgap!(fig.layout, 11)
    colgap!(fig.layout, 22)
    return fig
end

function main()
    length(ARGS) <= 1 || error("Expected at most one output path")
    output = if isempty(ARGS)
        joinpath(ROOT, "docs", "figures", "public_tg_method_checks.svg")
    else
        abspath(only(ARGS))
    end
    audit(name) = TOML.parsefile(joinpath(ROOT, "docs", "audits", "public_tg_$(name).toml"))
    save(
        output,
        method_figure(audit("friedman"), audit("vyazovkin"), audit("ramp_approximations")),
    )
    return println("Wrote $output")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
