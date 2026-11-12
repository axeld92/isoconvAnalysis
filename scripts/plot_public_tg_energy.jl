using CairoMakie
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const METHODS = (
    ("friedman", "Friedman"),
    ("advanced_vyazovkin", "Adv. Vyazovkin"),
    ("kas", "KAS"),
    ("fwo", "FWO"),
    ("starink", "Starink"),
)
const COLORS = ("#27638A", "#B64B25", "#367D62", "#5F5089", "#967727")
const MARKERS = (:circle, :rect, :diamond, :utriangle, :cross)

energy(profile) = [point["energy_kJ_per_mol"] for point in profile["points"]]

function selected_profile(profiles, method, case, protocol)
    return only(
        filter(profiles) do profile
            settings = profile["settings"]
            return profile["case"] == case &&
                   profile["method"] == method &&
                   get(settings, "derivative", protocol["primary_derivative"]) ==
                   protocol["primary_derivative"] &&
                   get(settings, "delta_alpha", protocol["primary_delta_alpha"]) ==
                   protocol["primary_delta_alpha"]
        end,
    )
end

function draw_series!(axis, alpha, values, index, label)
    return scatterlines!(
        axis,
        alpha,
        values;
        color=COLORS[index],
        marker=MARKERS[index],
        markersize=7,
        linewidth=1.6,
        label,
    )
end

function energy_figure(evidence, protocol)
    alpha = evidence["alpha_grid"]
    fig = Figure(; size=(1200, 810), fontsize=15, backgroundcolor=:white)
    constant = only(
        filter(c -> c["definition"]["id"] == "constant_affine", evidence["controls"])
    )
    variable = only(
        filter(c -> c["definition"]["id"] == "variable_affine", evidence["controls"])
    )
    constant_axis = Axis(
        fig[1, 1];
        title="Constant-energy affine control",
        xlabel="Conversion, α",
        ylabel="Estimate minus truth (kJ/mol)",
        xticks=0.1:0.2:0.9,
    )
    variable_axis = Axis(
        fig[1, 2];
        title="Conversion-dependent energy control",
        xlabel="Conversion, α",
        ylabel="Estimate minus pointwise truth (kJ/mol)",
        xticks=0.1:0.2:0.9,
    )
    hlines!(constant_axis, [0.0]; color="#888888", linewidth=1)
    hlines!(variable_axis, [0.0]; color="#888888", linewidth=1)
    for (i, (method, label)) in enumerate(METHODS)
        profile = selected_profile(
            constant["profiles"], method, "constant_affine", protocol
        )
        draw_series!(
            constant_axis,
            alpha,
            energy(profile) .- profile["truth_at_target_kJ_per_mol"],
            i,
            label,
        )
    end
    for (i, profile) in enumerate(variable["profiles"])
        label = if profile["method"] == "friedman"
            "Friedman"
        else
            "Vyazovkin Δα=$(profile["settings"]["delta_alpha"])"
        end
        draw_series!(
            variable_axis,
            alpha,
            energy(profile) .- profile["truth_at_target_kJ_per_mol"],
            i,
            label,
        )
    end
    Legend(
        fig[2, 1],
        constant_axis;
        orientation=:horizontal,
        nbanks=2,
        framevisible=false,
        labelsize=13,
        tellwidth=false,
    )
    Legend(
        fig[2, 2],
        variable_axis;
        orientation=:horizontal,
        nbanks=2,
        framevisible=false,
        labelsize=13,
        tellwidth=false,
    )

    full = only(filter(c -> c["id"] == "all_five", evidence["polypropylene_cohorts"]))
    omitted = only(filter(c -> c["id"] == "omit_10", evidence["polypropylene_cohorts"]))
    case = evidence["comparison_case"]
    pp_axis = Axis(
        fig[3, 1];
        title="Polypropylene, all five runs",
        xlabel="Conversion, α",
        ylabel="Apparent energy (kJ/mol)",
        xticks=0.1:0.2:0.9,
    )
    omit_axis = Axis(
        fig[3, 2];
        title="Effect of omitting the 10 K/min run",
        xlabel="Conversion, α",
        ylabel="Four-run minus five-run energy (kJ/mol)",
        xticks=0.1:0.2:0.9,
    )
    for (i, (method, _)) in enumerate(METHODS[1:2])
        scenarios = unique(p["case"] for p in full["profiles"] if p["method"] == method)
        profiles = [
            energy(selected_profile(full["profiles"], method, scenario, protocol)) for
            scenario in scenarios
        ]
        span = [extrema(values) for values in zip(profiles...)]
        band!(pp_axis, alpha, first.(span), last.(span); color=(COLORS[i], 0.15))
    end
    for (i, (method, label)) in enumerate(METHODS)
        primary = energy(selected_profile(full["profiles"], method, case, protocol))
        reduced = energy(selected_profile(omitted["profiles"], method, case, protocol))
        draw_series!(pp_axis, alpha, primary, i, label)
        draw_series!(omit_axis, alpha, reduced .- primary, i, label)
    end
    rowsize!(fig.layout, 2, Fixed(48))
    rowgap!(fig.layout, 12)
    colgap!(fig.layout, 25)
    return fig
end

function main()
    length(ARGS) <= 1 || error("Expected at most one output path")
    path = if isempty(ARGS)
        joinpath(ROOT, "docs", "figures", "public_tg_energy.svg")
    else
        abspath(only(ARGS))
    end
    evidence = TOML.parsefile(joinpath(ROOT, "docs", "audits", "public_tg_energy.toml"))
    protocol = TOML.parsefile(joinpath(ROOT, evidence["protocol"]))
    save(path, energy_figure(evidence, protocol))
    return println("Wrote $path")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
