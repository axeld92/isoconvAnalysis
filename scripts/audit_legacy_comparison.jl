module LegacyComparisonAudit

using Dates
using IsoconversionalAnalysis
using Printf
using SHA
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_REPORT = joinpath(PROJECT_ROOT, "docs", "legacy_comparison_audit.md")
const DEFAULT_EVIDENCE = joinpath(
    PROJECT_ROOT, "docs", "audits", "m9_legacy_comparison.toml"
)

const LEGACY_ALIAS_PATH = "code/2021mar.mat"
const LEGACY_CORRECTED_PATH = "code/2021marcorrected.mat"

const LEGACY_MAT_ARTIFACTS = [
    (
        path="code/2021corridas/2021corridas.mat",
        tier="canonical_raw",
        comparison_use="exact named-variable comparison",
    ),
    (
        path=LEGACY_ALIAS_PATH,
        tier="raw_alias",
        comparison_use="exact named-variable comparison",
    ),
    (
        path=LEGACY_CORRECTED_PATH,
        tier="cropped_raw_reference",
        comparison_use="exact contiguous-subset comparison",
    ),
    (
        path="code/isotherm.mat",
        tier="held_out_raw",
        comparison_use="Julia validation source; no separate legacy output oracle",
    ),
    (
        path="code/corridasjulio.mat",
        tier="derived_workspace",
        comparison_use="schema and workflow context only",
    ),
    (
        path="code/julio2.mat",
        tier="derived_workspace",
        comparison_use="schema and workflow context only",
    ),
    (
        path="code/000_nuevo.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/025_nuevo.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/050_nuevo.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/075_nuevo.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/100_nuevo.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/a00neum.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/a25neum.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/a50neum.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/a75neum.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/a100neum.mat",
        tier="saved_result_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/analisisarticuloneumatics.mat",
        tier="publication_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
    (
        path="code/articulo/articuloneumaticos.mat",
        tier="publication_workspace",
        comparison_use="snapshot context only; producing command is not recorded",
    ),
]

const METHOD_COMPARISONS = [
    (
        milestone="M2",
        area="Raw data ingestion",
        legacy="Direct MAT variables and St01-St20 aliases",
        julia="Checksum-verified, explicitly named MAT variables with preserved invalid rows",
        status="exact",
        conclusion="All 20 configured dynamic arrays are exactly equal to their aliases, including NaNs.",
    ),
    (
        milestone="M3",
        area="Preprocessing",
        legacy="Row endpoints, Savitzky-Golay smoothing, diff_alpha, and two low-pass calls",
        julia="Explicit 150-700 C references, local-polynomial smoothing, diagnostics, and eight sensitivity profiles",
        status="intentional_deviation",
        conclusion="The scientific conversion-and-derivative objective is preserved; numerical parity is neither available nor accepted.",
    ),
    (
        milestone="M4",
        area="Linear isoconversional methods",
        legacy="Four-run Friedman, KAS, FWO, and Starink scripts",
        julia="Named multi-experiment methods with unit checks, first-crossing interpolation, and two-sided slope intervals",
        status="corrected_equations",
        conclusion="Method families are preserved, while known Friedman and temperature-handling defects are corrected.",
    ),
    (
        milestone="M4",
        area="Advanced Vyazovkin method",
        legacy="Optimization of distance from a hard-coded objective value of 12",
        julia="Bounded minimization of the true pairwise objective with Fisher intervals",
        status="corrected_equations",
        conclusion="The method concept is preserved; the corrected objective is expected to produce different values.",
    ),
    (
        milestone="M5",
        area="Fraser-Suzuki deconvolution",
        legacy="Global/workspace fits, complex-to-zero masking, and sequentially tightened bounds",
        julia="Stable real-valued peak definition, independent deterministic multistart fits, and structural model selection",
        status="intentional_deviation",
        conclusion="Curve decomposition is preserved, but order dependence and invalid-domain masking are removed.",
    ),
    (
        milestone="M6",
        area="Kinetic triplets and reaction models",
        legacy="Normal-equation regressions with duplicated or confounded intercept terms",
        julia="Single-intercept QR fits, a named model registry, diagnostics, and held-out-rate gates",
        status="intentional_deviation",
        conclusion="The candidate equations are preserved as empirical models; legacy coefficients are not valid parameter oracles.",
    ),
    (
        milestone="M7",
        area="Dynamic and ramp-and-hold prediction",
        legacy="Workspace-dependent simulations with hard-coded runs and no frozen predictive gate",
        julia="Bounded integration, calibration-only parameters, and 15 held-out ramp-and-hold programs",
        status="not_numerically_comparable",
        conclusion="Both evaluated Julia formulations are rejected; matching legacy plots would not establish prediction quality.",
    ),
    (
        milestone="M8",
        area="Mixture additivity",
        legacy="Observed-minus-additive conversion/rate curves after legacy preprocessing",
        julia="Matched-rate measured curves, three channels, signed integrals, and preprocessing envelopes",
        status="concept_preserved",
        conclusion="The sign convention and additive concept agree, but corrected preprocessing and expanded diagnostics preclude pointwise parity.",
    ),
]

const INTENTIONAL_DEVIATIONS = [
    (
        id="LD01",
        area="Data schema",
        legacy_evidence="Raw arrays use mass in column 4; derived five-column workspaces use column 3.",
        julia_treatment="Select schemas by catalog source and named columns, never by column count alone.",
        reason="Prevents silent interpretation of DSC or derivative values as mass.",
    ),
    (
        id="LD02",
        area="Corrected workspace",
        legacy_evidence="Every 2021marcorrected array is an exact contiguous copy of raw columns 1-4 beginning near 110 C.",
        julia_treatment="Classify it as a cropped raw reference, not processed ground truth.",
        reason="The artifact contains no normalized conversion or smoothed derivative to reproduce.",
    ),
    (
        id="LD03",
        area="Preprocessing",
        legacy_evidence="cleandata.m normalizes with retained endpoint rows and applies undocumented smoothing and low-pass choices.",
        julia_treatment="Estimate exact-temperature reference masses robustly and expose smoothing and endpoint sensitivity.",
        reason="Makes the conversion coordinate and derivative bias explicit and testable.",
    ),
    (
        id="LD04",
        area="Friedman regression",
        legacy_evidence="friedman.m has incorrect heating-rate parentheses and interpolates Temp_4 instead of da_4.",
        julia_treatment="Use measured heating rate and the conversion derivative for every experiment.",
        reason="Never preserve known algebraic or variable-selection defects for parity.",
    ),
    (
        id="LD05",
        area="Isoconversional generality",
        legacy_evidence="Several scripts accept exactly four runs and mix Celsius-to-kelvin conversion conventions.",
        julia_treatment="Accept arbitrary valid experiment groups and enforce kelvin at the typed boundary.",
        reason="Removes fixed-design and double-conversion failure modes.",
    ),
    (
        id="LD06",
        area="Advanced Vyazovkin objective",
        legacy_evidence="Vyazovkin.m minimizes abs(phiofe2(...)-12), where 12 assumes four equal pair terms.",
        julia_treatment="Minimize the true pairwise objective with bounded search and method-specific intervals.",
        reason="The optimum must not depend on a hard-coded experiment-count target.",
    ),
    (
        id="LD07",
        area="Fraser-Suzuki domain",
        legacy_evidence="deconvolve.m replaces complex-valued evaluations with zero.",
        julia_treatment="Use a stable Gaussian limit and a continuous real-valued zero outside the logarithm domain.",
        reason="Avoids silently accepting invalid intermediate arithmetic.",
    ),
    (
        id="LD08",
        area="Deconvolution reproducibility",
        legacy_evidence="deconvolve.m hard-codes St05; deconvolution4.m tightens bounds from each preceding fit.",
        julia_treatment="Fit each named experiment independently with deterministic multistart optimization.",
        reason="Removes hidden workspace inputs and experiment-order dependence.",
    ),
    (
        id="LD09",
        area="Peak accounting and interpretation",
        legacy_evidence="deconvolve4.m omits peak 4 from denominators and assigns its fraction from peak 3.",
        julia_treatment="Integrate every selected peak and keep peaks ordered but chemically unnamed.",
        reason="Corrects area fractions without inventing component identity.",
    ),
    (
        id="LD10",
        area="Kinetic regression",
        legacy_evidence="kineticanalysis.m duplicates the constant column and solves least squares through normal equations.",
        julia_treatment="Fit one ln(A/min^-1) intercept using QR and retain rank, conditioning, and predictive diagnostics.",
        reason="Restores parameter identifiability and numerical stability.",
    ),
    (
        id="LD11",
        area="Reaction-model intervals",
        legacy_evidence="sblin.m retains amplitude/prefactor confounding; sbrlin.m reports observation-sized prediction widths as coefficient intervals.",
        julia_treatment="Separate coefficient intervals, residual diagnostics, model criteria, and cross-rate prediction.",
        reason="Avoids overstating parameter certainty from densely sampled derivative points.",
    ),
    (
        id="LD12",
        area="Compensation and simulation",
        legacy_evidence="compeffects.m includes a wrong Ginstling-Brounshtein expression; the simulation path uses globals and hard-coded St16.",
        julia_treatment="Keep compensation diagnostic and simulate only frozen, provenance-carrying, bounded formulations.",
        reason="Prevents circular parameter substitution and nonphysical hidden-state propagation.",
    ),
    (
        id="LD13",
        area="Mixture interaction",
        legacy_evidence="synergy.m uses observed-minus-additive differences after legacy preprocessing.",
        julia_treatment="Preserve the sign convention while reporting mass, conversion, rate, signed areas, and preprocessing scenarios.",
        reason="Retains the scientific comparison without claiming a named chemical synergy mechanism.",
    ),
]

portable_path(path::AbstractString) = replace(normpath(path), '\\' => '/')

function absolute_path(root::AbstractString, relative::AbstractString)
    return normpath(joinpath(root, relative))
end

function sha256_file(path::AbstractString)
    return open(path, "r") do io
        return bytes2hex(sha256(io))
    end
end

function require_condition(condition::Bool, message::AbstractString)
    condition || throw(ArgumentError("legacy-comparison audit failed: $message"))
    return nothing
end

function exact_contiguous_bounds(raw::AbstractMatrix, cropped::AbstractMatrix)
    size(cropped, 2) <= size(raw, 2) || return nothing
    size(cropped, 1) <= size(raw, 1) || return nothing
    column_count = size(cropped, 2)
    candidate = findfirst(axes(raw, 1)) do row
        return all(
            column -> isequal(raw[row, column], cropped[1, column]), axes(cropped, 2)
        )
    end
    isnothing(candidate) && return nothing
    final_row = candidate + size(cropped, 1) - 1
    final_row <= size(raw, 1) || return nothing
    raw_view = @view raw[candidate:final_row, 1:column_count]
    all(isequal.(raw_view, cropped)) || return nothing
    return (candidate, final_row)
end

function collect_experiment_evidence(root::AbstractString=PROJECT_ROOT)
    catalog = load_dataset_catalog(joinpath(root, "config", "datasets.toml"))
    dynamic_source = only(filter(source -> source.id == "dynamic_2021", catalog.sources))
    alias_path = absolute_path(root, LEGACY_ALIAS_PATH)
    corrected_path = absolute_path(root, LEGACY_CORRECTED_PATH)
    rows = NamedTuple[]
    for experiment in dynamic_source.experiments
        isnothing(experiment.legacy_alias) && continue
        canonical = load_mat_variable(
            dynamic_source.source_file, experiment.source_variable
        )
        alias = load_mat_variable(alias_path, experiment.legacy_alias)
        corrected = load_mat_variable(corrected_path, experiment.legacy_alias)
        bounds = exact_contiguous_bounds(canonical, corrected)
        push!(
            rows,
            (
                experiment_id=experiment.id,
                canonical_variable=experiment.source_variable,
                legacy_alias=experiment.legacy_alias,
                canonical_rows=size(canonical, 1),
                canonical_columns=size(canonical, 2),
                alias_rows=size(alias, 1),
                alias_columns=size(alias, 2),
                alias_exact=isequal(canonical, alias),
                corrected_rows=size(corrected, 1),
                corrected_columns=size(corrected, 2),
                corrected_exact_contiguous_subset=(!isnothing(bounds)),
                corrected_start_row=isnothing(bounds) ? 0 : bounds[1],
                corrected_end_row=isnothing(bounds) ? 0 : bounds[2],
                corrected_start_temperature_celsius=corrected[1, 1],
                corrected_end_temperature_celsius=corrected[end, 1],
            ),
        )
    end
    sort!(rows; by=row -> row.experiment_id)
    return rows
end

function collect_artifact_evidence(root::AbstractString=PROJECT_ROOT)
    artifacts = NamedTuple[]
    for spec in LEGACY_MAT_ARTIFACTS
        path = absolute_path(root, spec.path)
        require_condition(isfile(path), "missing $(spec.path)")
        push!(
            artifacts,
            (
                path=spec.path,
                tier=spec.tier,
                comparison_use=spec.comparison_use,
                bytes=filesize(path),
                sha256=sha256_file(path),
            ),
        )
    end
    sort!(artifacts; by=artifact -> artifact.path)
    return artifacts
end

function validate_evidence(experiments, artifacts)
    require_condition(length(experiments) == 20, "expected 20 dynamic comparisons")
    require_condition(
        all(row -> row.alias_exact, experiments), "at least one raw alias differs"
    )
    require_condition(
        all(row -> row.corrected_exact_contiguous_subset, experiments),
        "at least one corrected array is not an exact raw subset",
    )
    require_condition(
        all(row -> row.corrected_columns == 4, experiments),
        "corrected arrays do not all have four columns",
    )
    require_condition(
        length(unique(row.id for row in INTENTIONAL_DEVIATIONS)) ==
        length(INTENTIONAL_DEVIATIONS),
        "intentional-deviation ids are not unique",
    )
    require_condition(
        length(artifacts) == length(LEGACY_MAT_ARTIFACTS),
        "artifact inventory is incomplete",
    )
    return nothing
end

function render_report(experiments, artifacts; generated_on::Date=Dates.today())
    validate_evidence(experiments, artifacts)
    alias_count = count(row -> row.alias_exact, experiments)
    crop_count = count(row -> row.corrected_exact_contiguous_subset, experiments)
    minimum_start = minimum(row.corrected_start_temperature_celsius for row in experiments)
    maximum_start = maximum(row.corrected_start_temperature_celsius for row in experiments)
    minimum_end = minimum(row.corrected_end_temperature_celsius for row in experiments)
    maximum_end = maximum(row.corrected_end_temperature_celsius for row in experiments)
    saved_workspace_count = count(
        artifact -> artifact.tier in ("saved_result_workspace", "publication_workspace"),
        artifacts,
    )

    io = IOBuffer()
    println(io, "# M9 corrected Julia-versus-legacy comparison")
    println(io)
    println(
        io,
        "Generated by `scripts/audit_legacy_comparison.jl` on $(generated_on). The audit " *
        "compares preserved legacy evidence with the corrected Julia workflow; it does not " *
        "treat every historical workspace value as ground truth.",
    )
    println(io)
    println(io, "## Decision")
    println(io)
    println(
        io,
        "Exact numerical parity is required and achieved at the raw dynamic import boundary: " *
        "$alias_count/20 canonical variables are exactly equal to `St01`-`St20`, including " *
        "the trailing NaN rows. Beyond that boundary, the Julia implementation preserves " *
        "scientific method families while intentionally correcting identified defects, " *
        "making assumptions explicit, and applying synthetic or held-out validation.",
    )
    println(io)
    println(
        io,
        "The file named `2021marcorrected.mat` is not a numerical preprocessing oracle. All " *
        "$crop_count/20 arrays are exact contiguous copies of raw columns 1-4, with start " *
        @sprintf(
            "temperatures %.3f-%.3f C and end temperatures %.3f-%.3f C. It contains no stored normalized conversion, smoothed derivative, fit, or uncertainty result.",
            minimum_start,
            maximum_start,
            minimum_end,
            maximum_end,
        ),
    )
    println(io)
    println(
        io,
        "The $saved_workspace_count saved result/publication workspaces have no recorded " *
        "producing command, environment, script revision, or unambiguous variable lineage in " *
        "the archive. They remain qualitative workflow evidence and are excluded from release " *
        "acceptance thresholds.",
    )
    println(io)
    println(io, "## Exact raw and cropped-array evidence")
    println(io)
    println(
        io,
        "| Experiment | Canonical variable | Alias | Raw shape | Exact alias | Corrected rows | Raw row interval | Corrected temperature range (C) |",
    )
    println(io, "|---|---|---|---:|:---:|---:|---:|---:|")
    for row in experiments
        raw_shape = "$(row.canonical_rows)x$(row.canonical_columns)"
        interval = "$(row.corrected_start_row):$(row.corrected_end_row)"
        temperatures = @sprintf(
            "%.3f-%.3f",
            row.corrected_start_temperature_celsius,
            row.corrected_end_temperature_celsius,
        )
        println(
            io,
            "| `$(row.experiment_id)` | `$(row.canonical_variable)` | " *
            "`$(row.legacy_alias)` | $raw_shape | $(row.alias_exact) | " *
            "$(row.corrected_rows) | $interval | $temperatures |",
        )
    end
    println(io)
    println(io, "## Method-level comparison")
    println(io)
    println(
        io, "| Milestone | Area | Legacy behavior | Julia behavior | Status | Conclusion |"
    )
    println(io, "|---|---|---|---|---|---|")
    for row in METHOD_COMPARISONS
        println(
            io,
            "| $(row.milestone) | $(row.area) | $(row.legacy) | $(row.julia) | " *
            "`$(row.status)` | $(row.conclusion) |",
        )
    end
    println(io)
    println(io, "## Intentional deviations")
    println(io)
    println(io, "| ID | Area | Legacy evidence | Julia treatment | Reason |")
    println(io, "|---|---|---|---|---|")
    for row in INTENTIONAL_DEVIATIONS
        println(
            io,
            "| `$(row.id)` | $(row.area) | $(row.legacy_evidence) | " *
            "$(row.julia_treatment) | $(row.reason) |",
        )
    end
    println(io)
    println(io, "## Evidence tiers and provenance")
    println(io)
    println(
        io,
        "The machine record at `docs/audits/m9_legacy_comparison.toml` stores a SHA-256 hash " *
        "for every MAT artifact used in this comparison. The tiers have " *
        "different evidentiary force:",
    )
    println(io)
    println(
        io,
        "- `canonical_raw` and `raw_alias`: exact numerical comparison is authoritative.",
    )
    println(
        io,
        "- `cropped_raw_reference`: authoritative only for the exact raw subset it contains.",
    )
    println(
        io,
        "- `held_out_raw`: authoritative measured input, reserved for prediction validation.",
    )
    println(
        io,
        "- `derived_workspace`: schema and workflow context; transformation lineage is incomplete.",
    )
    println(
        io,
        "- `saved_result_workspace` and `publication_workspace`: snapshot context only, not numerical acceptance oracles.",
    )
    println(
        io,
        "The method comparison preserves findings from the historical script review. " *
        "The original MATLAB source is stored separately.",
    )
    println(io)
    println(io, "## Release consequence")
    println(io)
    println(
        io,
        "The M9 legacy comparison and intentional-deviation record are complete. A future " *
        "difference in raw alias equality is a regression. A difference in corrected " *
        "preprocessing, fitted activation energies, peak parameters, reaction-model " *
        "coefficients, or simulations is not by itself a regression; those results are " *
        "governed by the Julia equations, diagnostics, frozen configurations, synthetic " *
        "tests, and predictive gates documented in M3-M8.",
    )
    println(io)
    println(
        io,
        "Subsequent release work is tracked in PROJECT_STATUS.md; this audit remains the " *
        "frozen evidence for the legacy-comparison and intentional-deviation portion of M9.",
    )
    return String(take!(io))
end

function toml_quote(value::AbstractString)
    escaped = replace(value, '\\' => "\\\\", '"' => "\\\"")
    return "\"$escaped\""
end

function render_machine_evidence(
    root::AbstractString,
    experiments,
    artifacts,
    report_path::AbstractString;
    generated_on::Date=Dates.today(),
)
    validate_evidence(experiments, artifacts)
    report_hash = sha256_file(report_path)
    generator_hash = sha256_file(joinpath(root, "scripts", "audit_legacy_comparison.jl"))
    fingerprint_material = join(
        vcat(
            [report_hash, generator_hash],
            [artifact.sha256 for artifact in artifacts],
            [
                join(
                    (
                        row.experiment_id,
                        string(row.alias_exact),
                        string(row.corrected_exact_contiguous_subset),
                        string(row.corrected_start_row),
                        string(row.corrected_end_row),
                    ),
                    ":",
                ) for row in experiments
            ],
            [row.id * ":" * row.julia_treatment for row in INTENTIONAL_DEVIATIONS],
        ),
        "|",
    )
    analysis_fingerprint = bytes2hex(sha256(fingerprint_material))

    io = IOBuffer()
    println(io, "schema_version = 1")
    println(io, "generated_on = ", toml_quote(string(generated_on)))
    println(io, "generator = ", toml_quote("scripts/audit_legacy_comparison.jl"))
    println(io, "generator_sha256 = ", toml_quote(generator_hash))
    println(io, "report = ", toml_quote(portable_path(relpath(report_path, root))))
    println(io, "report_sha256 = ", toml_quote(report_hash))
    println(io, "analysis_fingerprint = ", toml_quote(analysis_fingerprint))
    println(io, "comparison_scope = ", toml_quote("corrected_julia_versus_legacy"))
    println(io, "numerical_parity_scope = ", toml_quote("raw_dynamic_import_only"))
    println(io, "legacy_result_workspaces_are_numerical_oracles = false")
    println(io, "dynamic_experiment_count = ", length(experiments))
    println(io, "exact_raw_alias_count = ", count(row -> row.alias_exact, experiments))
    println(
        io,
        "exact_corrected_contiguous_subset_count = ",
        count(row -> row.corrected_exact_contiguous_subset, experiments),
    )
    println(io, "method_comparison_count = ", length(METHOD_COMPARISONS))
    println(io, "intentional_deviation_count = ", length(INTENTIONAL_DEVIATIONS))
    println(io, "artifact_count = ", length(artifacts))

    for row in experiments
        println(io)
        println(io, "[[experiment]]")
        println(io, "id = ", toml_quote(row.experiment_id))
        println(io, "canonical_variable = ", toml_quote(row.canonical_variable))
        println(io, "legacy_alias = ", toml_quote(row.legacy_alias))
        println(io, "canonical_rows = ", row.canonical_rows)
        println(io, "canonical_columns = ", row.canonical_columns)
        println(io, "alias_rows = ", row.alias_rows)
        println(io, "alias_columns = ", row.alias_columns)
        println(io, "alias_exact_including_nan = ", row.alias_exact)
        println(io, "corrected_rows = ", row.corrected_rows)
        println(io, "corrected_columns = ", row.corrected_columns)
        println(
            io,
            "corrected_exact_contiguous_subset = ",
            row.corrected_exact_contiguous_subset,
        )
        println(io, "corrected_start_row_1_based = ", row.corrected_start_row)
        println(io, "corrected_end_row_1_based = ", row.corrected_end_row)
        println(
            io,
            "corrected_start_temperature_celsius = ",
            repr(row.corrected_start_temperature_celsius),
        )
        println(
            io,
            "corrected_end_temperature_celsius = ",
            repr(row.corrected_end_temperature_celsius),
        )
    end

    for row in METHOD_COMPARISONS
        println(io)
        println(io, "[[method_comparison]]")
        println(io, "milestone = ", toml_quote(row.milestone))
        println(io, "area = ", toml_quote(row.area))
        println(io, "legacy = ", toml_quote(row.legacy))
        println(io, "julia = ", toml_quote(row.julia))
        println(io, "status = ", toml_quote(row.status))
        println(io, "conclusion = ", toml_quote(row.conclusion))
    end

    for row in INTENTIONAL_DEVIATIONS
        println(io)
        println(io, "[[intentional_deviation]]")
        println(io, "id = ", toml_quote(row.id))
        println(io, "area = ", toml_quote(row.area))
        println(io, "legacy_evidence = ", toml_quote(row.legacy_evidence))
        println(io, "julia_treatment = ", toml_quote(row.julia_treatment))
        println(io, "reason = ", toml_quote(row.reason))
    end

    for artifact in artifacts
        println(io)
        println(io, "[[artifact]]")
        println(io, "path = ", toml_quote(artifact.path))
        println(io, "tier = ", toml_quote(artifact.tier))
        println(io, "comparison_use = ", toml_quote(artifact.comparison_use))
        println(io, "bytes = ", artifact.bytes)
        println(io, "sha256 = ", toml_quote(artifact.sha256))
    end
    return String(take!(io))
end

function write_text(path::AbstractString, content::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        return print(io, content)
    end
    return path
end

function main(args=ARGS)
    isempty(args) || throw(ArgumentError("this audit does not accept arguments"))
    experiments = collect_experiment_evidence(PROJECT_ROOT)
    artifacts = collect_artifact_evidence(PROJECT_ROOT)
    generated_on = Dates.today()
    report = render_report(experiments, artifacts; generated_on=generated_on)
    write_text(DEFAULT_REPORT, report)
    machine = render_machine_evidence(
        PROJECT_ROOT, experiments, artifacts, DEFAULT_REPORT; generated_on=generated_on
    )
    write_text(DEFAULT_EVIDENCE, machine)
    TOML.parsefile(DEFAULT_EVIDENCE)
    println("Wrote $(portable_path(relpath(DEFAULT_REPORT, PROJECT_ROOT)))")
    println("Wrote $(portable_path(relpath(DEFAULT_EVIDENCE, PROJECT_ROOT)))")
    return nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    LegacyComparisonAudit.main(ARGS)
end
