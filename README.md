# IsoconversionalAnalysis.jl

A from-scratch Julia implementation of the thermogravimetric and isoconversional-analysis
ideas developed in the original master thesis. The legacy MATLAB directory is retained as
read-only evidence; undocumented legacy behavior is not treated as a specification.

The concise current snapshot and next actions are in
[`PROJECT_STATUS.md`](PROJECT_STATUS.md). The complete scientific scope, decisions, and
milestone history live in [`MASTERPLAN.md`](MASTERPLAN.md).
M0 documents the source data and methods, M1 establishes the reproducible Julia foundation,
M2 imports and audits the source data, M3 provides the validated preprocessing pipeline, M4
implements and audits five isoconversional methods, and M5 provides constrained Fraser–Suzuki
deconvolution with peak-count and identifiability diagnostics. M6 adds a named differential
reaction-model registry, identifiable kinetic-triplet regression, model comparison,
heating-rate-level validation, and a guarded compensation-effect diagnostic.

## Quick start

Julia 1.11 is required.

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. -e 'using IsoconversionalAnalysis; println("loaded")'
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

Load the validated default configuration with:

```julia
using IsoconversionalAnalysis

config = load_config("config/analysis_defaults.toml")
```

Load the complete M2 experiment catalog with:

```julia
catalog = load_dataset_catalog("config/datasets.toml")
experiments = load_experiments(catalog)
```

The importer verifies source checksums and loads only explicitly named MAT variables. Run
`julia --startup-file=no --project=. scripts/audit_legacy_data.jl` to reproduce the raw-data
inventory.

Preprocess the calibration experiments with the accepted M3 profile:

```julia
processed = preprocess.(experiments, Ref(config))
```

Here `experiments` should contain the calibration runs loaded above. Regenerate the M3
sensitivity report and static QC figures with:

```sh
julia --startup-file=no --project=docs scripts/audit_preprocessing.jl
```

The accepted settings, alternatives, warnings, and real-data results are documented in
[`docs/preprocessing_audit.md`](docs/preprocessing_audit.md).

Group processed calibration runs by composition and estimate activation-energy profiles with
the shared M4 API:

```julia
group = filter(experiment -> experiment.composition.waste_tire == 0.5, processed)
starink = analyze(:starink, group, config)
advanced = analyze(:advanced_vyazovkin, group, config)
```

Supported method symbols are `:friedman`, `:kas`, `:fwo`, `:starink`, and
`:advanced_vyazovkin`. Regenerate the complete real-data and endpoint-sensitivity audit with:

```sh
julia --startup-file=no --project=docs scripts/audit_isoconversional.jl
```

The conventions, diagnostics, figures, and complete primary results are documented in
[`docs/isoconversional_audit.md`](docs/isoconversional_audit.md).

Compare two-, three-, and four-peak Fraser–Suzuki decompositions for one processed run with:

```julia
comparison = compare_peak_counts(first(group), config)
fit = only(filter(result -> result.peak_count == comparison.selected_peak_count,
                  comparison.results))
```

The raw information-criterion minimum is retained separately from the structurally accepted
selection. Regenerate all 20 independent comparisons and the five evaluated joint multi-rate
fits with:

```sh
julia --startup-file=no --project=docs scripts/audit_deconvolution.jl
```

The interpretation limits, diagnostics, figures, and complete results are documented in
[`docs/deconvolution_audit.md`](docs/deconvolution_audit.md).

Compare all configured M6 reaction models jointly across the four dynamic rates for a
composition with:

```julia
comparison = compare_reaction_models(group, config)
fit = only(filter(result -> result.model.name == comparison.selected_model,
                  comparison.results))
```

Regenerate the five composition comparisons, all leave-one-heating-rate-out validations, and
the guarded compensation diagnostics with:

```sh
julia --startup-file=no --project=docs scripts/audit_reaction_models.jl
```

The complete results and predictive decision are documented in
[`docs/reaction_model_audit.md`](docs/reaction_model_audit.md). The real-data gate rejects a
constant kinetic triplet as the default predictive description; M7 will retain it as a
benchmark while evaluating conversion-dependent parameters.

Development, formatting, and documentation commands are listed in
[`docs/src/development.md`](docs/src/development.md). Dependency rationale and licenses are
recorded in [`docs/dependencies.md`](docs/dependencies.md).
