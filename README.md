# IsoconversionalAnalysis.jl

IsoconversionalAnalysis.jl is a Julia toolkit for analyzing thermogravimetric data.
It estimates how apparent activation energy changes with conversion, compares
candidate kinetic models, and records the assumptions behind the results.

## Capabilities

- Import segment-aware TG data from declared columns and units with source checksum verification.
- Preprocess mass-loss curves with configurable reference masses, conversion
  windows, smoothing, and derivative settings.
- Estimate apparent activation energy with Friedman, KAS, Flynn–Wall–Ozawa,
  Starink, and advanced Vyazovkin methods.
- Fit overlapping Fraser–Suzuki peaks and compare candidate reaction models and
  kinetic triplets. Inspect diagnostics and compare fits across heating rates.

## Quick start

From a clone of this repository, use Julia 1.11 to install the dependencies.

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
```

Start Julia with `julia --project=.`. The following example estimates an apparent
activation-energy profile for one sample group in the bundled thermogravimetry
dataset.

```julia
using IsoconversionalAnalysis

config = load_config("config/analysis_defaults.toml")
catalog = load_dataset_catalog("config/datasets.toml")
experiments = load_experiments(catalog; roles=[:calibration])
group = filter(experiment -> experiment.composition.waste_tire == 0.5, experiments)
curves = preprocess.(group, Ref(config))
profile = analyze(:advanced_vyazovkin, curves, config)
```

Run the tests from the repository root.

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

- [Data ingestion and provenance](docs/src/data_ingestion.md)
- [Preprocessing and sensitivity checks](docs/preprocessing_audit.md)
- [Activation-energy methods and diagnostics](docs/src/isoconversional.md)
- [Peak fitting](docs/src/deconvolution.md) and [reaction-model comparison](docs/src/reaction_models.md)
- [Development guide](docs/src/development.md) and [data-sharing status](DATA_SHARING.md)

## Public TG benchmark

Reproduce known-energy recovery and the polypropylene activation-energy profiles
using the [public benchmark instructions](docs/public_tg_energy.md). The
[TG workflow guide](docs/src/thermal_conversion.md) explains how to check common
conversion crossings, temperature coverage, local rates, recorded-time intervals,
and ramp assumptions.

## Research direction

The project is developing DSC reaction-progress workflows, model fitting across
multiple experiments, uncertainty assessment, and prediction
on thermal programs excluded from fitting.
