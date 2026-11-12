using CairoMakie
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const GROUPS = (
    ("Constant_E", "Constant-energy synthetic TG"),
    ("Two_Steps", "Two-step synthetic TG"),
    ("Variable_E", "Variable-energy synthetic TG"),
    ("polypropylene", "Experimental polypropylene TG"),
)

function conversion_support_figure(evidence)
    fig = Figure(; size=(1100, 770), fontsize=16, backgroundcolor=:white)
    for (i, (material, title)) in enumerate(GROUPS)
        group = only(filter(g -> g["material_id"] == material, evidence["groups"]))
        case = only(filter(c -> c["id"] == evidence["comparison_case"], group["cases"]))
        points = filter(p -> p["common_upward"], case["points"])
        x = [p["alpha"] for p in points]
        axis = Axis(
            fig[cld(i, 2), isodd(i) ? 1 : 2];
            title="$title ($(length(group["run_ids"])) runs)",
            xlabel="Conversion, α",
            ylabel="Crossing temperature (°C)",
            xticks=0.0:0.2:1.0,
        )
        xlims!(axis, -0.025, 1.025)
        temperature(field) = [p[field] - 273.15 for p in points]
        bars = rangebars!(
            axis,
            x,
            temperature("temperature_min_K"),
            temperature("temperature_max_K");
            color="#356B8A",
            linewidth=2.5,
            whiskerwidth=8,
        )
        means = scatter!(
            axis, x, temperature("temperature_mean_K"); color="#A84123", markersize=8
        )
        if i == length(GROUPS)
            Legend(
                fig[3, 1:2],
                [bars, means],
                ["Min–max across runs", "Mean crossing temperature"];
                orientation=:horizontal,
                framevisible=false,
                labelsize=14,
            )
        end
    end
    rowsize!(fig.layout, 3, Fixed(35))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 22)
    return fig
end

function main()
    length(ARGS) <= 1 || error("Expected at most one output path")
    output = if isempty(ARGS)
        joinpath(ROOT, "docs", "figures", "public_tg_conversion_support.svg")
    else
        abspath(only(ARGS))
    end
    path = joinpath(ROOT, "docs", "audits", "public_tg_conversion_support.toml")
    save(output, conversion_support_figure(TOML.parsefile(path)))
    return println("Wrote $output")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
