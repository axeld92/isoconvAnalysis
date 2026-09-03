# Data ingestion and audit

M2 treats data loading as a strict, read-only boundary. `config/datasets.toml` explicitly
maps 35 stable experiment ids to two fingerprinted MAT v5 files and named variables:

- 20 dynamic calibration runs from `code/2021corridas/2021corridas.mat`;
- 15 ramp-and-hold validation runs from `code/isotherm.mat`.

The corresponding NETZSCH text-export path is retained for every dynamic run. The MAT
snapshot is the M2 import boundary; the instrument exports remain the primary provenance
source described in `docs/data_dictionary.md`.

## Loading

```julia
using IsoconversionalAnalysis

catalog = load_dataset_catalog("config/datasets.toml")
all_experiments = load_experiments(catalog)
calibration_runs = load_experiments(catalog; roles=[:calibration])
held_out_runs = load_experiments(catalog; roles=[:validation])
one_run = load_experiment(catalog, "dynamic_wt050_rate10")
```

Every source checksum is verified before values are read. Only explicitly configured MAT
variables are loaded. Missing variables, checksum mismatches, wrong matrix types, unexpected
column counts, duplicate catalog ids, and invalid metadata raise dedicated exceptions.

## Raw `Experiment` contract

`Experiment` stores temperature in kelvin, time in minutes, remaining mass in percent, the
available DSC/flow/sensitivity/segment channels, dry-mass composition, declared conditions,
source hash, and provenance metadata.

M2 does not crop or repair arrays. In particular, the trailing all-`NaN` row present in every
selected MATLAB matrix is retained and recorded in `import_warnings`. Segment values that are
absent are represented as `nothing`; invalid individual segment values are represented as
`missing`. M3 will make the explicit row-removal and interval-selection decisions.

Use `validate_experiment`, `valid_row_mask`, and `trailing_invalid_row_count` to inspect the
raw structural state.

## Reproducing the audit

```sh
julia --startup-file=no --project=. scripts/audit_legacy_data.jl
```

This regenerates:

- `docs/data_audit.md`, the human-readable audit;
- `docs/audits/m2_data_inventory.toml`, the machine-readable inventory.

The audit reports ranges, row validity, sampling intervals, ordering, segment values, and
measured heating-rate slopes. Dynamic slopes use 100–700 °C; ramp-and-hold slopes use recorded
segment 3. These are diagnostics and do not replace M3 preprocessing.
