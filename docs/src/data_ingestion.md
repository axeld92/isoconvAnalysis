# Data ingestion and audit

## Delimited TG and DSC observations

Use `read_thermal_delimited` to read CSV, TSV, or semicolon-separated exports with
declared column names and units. For example, the experimental polypropylene files in
the public benchmark cache use tabs despite their `.csv` suffixes:

```julia
using IsoconversionalAnalysis

mapping = ThermalColumnMapping(
    ThermalColumnSpec("Temperature (C)", :C),
    [ThermalSignalSpec(:tg, ThermalColumnSpec("Weight (mg)", :mg), :mass)];
    time=ThermalColumnSpec("Time (min)", :min),
)
run = read_thermal_delimited("path/to/polypropylene.csv", mapping;
    material_id="polypropylene", delimiter='\t')
T = temperature_K(run)
t = time_min(run)
mass_mg = signal_values(run, :tg)
```

Map DSC channels as `:heat_flow` for power (such as `:mW`) or `:specific_heat_flow`
for power per mass (such as `:mW_per_mg`). Set `exotherm_sign=-1` for an export that
records exotherms downward. The sign convention accompanies the unchanged signal.
Multiple channels can be retained in one run, including simultaneously recorded TG,
DSC, and auxiliary observations. For the TURBO exports, use `encoding=:cp1252`,
`header_prefix="##"`, `delimiter=';'`, `segment_column="Segment"`, and the literal
temperature header `Temp./›C` after removing the prefix. The public manifest lists the
channel mappings for the dynamic and isothermal file variants.

For a named header followed by a separate row of units, pass `unit_labels` with
the literal labels in source column order. The importer validates that row and
keeps the original line numbers for the numerical observations. The cast-PMMA
exports use `unit_labels=["[#]", "[°C]", "[s]", "[mg]", "[°C]"]`.

`run.segments` lists contiguous blocks in acquisition order. Use
`select_thermal_segment(run, ordinal)` after identifying the desired phase from the
recorded label, metadata, and temperature history. A repeated label can occur in several
blocks. `thermal_segment_diagnostics(run)` reports temperature range, recorded duration,
temperature steps, and the descriptive T–t slope for each block. With no segment column,
the run contains one unlabelled block. With no time column, duration and slope are absent.
Recorded times must increase within a block. The importer retains clock resets across
blocks and temperature reversals within them.

The run retains native values and units. `temperature_K`, `time_min`, and `signal_values`
return independent arrays with optional unit conversion for signals. Converting mass to
relative mass requires a reference mass. Converting heat flow to reaction progress
requires an observation model for baseline and enthalpy. Specify `expected_sha256` to
verify source bytes before parsing. `run.provenance` identifies the source hash, physical
line numbers, and unmapped headers. Segment selections carry the source lines for their
subset. Malformed mapped data raises `DataImportError` with source context.

## Historical MAT catalog

M2 treats data loading as a strict, read-only boundary. `config/datasets.toml` explicitly
maps 35 stable experiment ids to two fingerprinted MAT v5 files and named variables:

- 20 dynamic calibration runs from `code/2021corridas/2021corridas.mat`;
- 15 ramp-and-hold validation runs from `code/isotherm.mat`.

The corresponding NETZSCH text-export path is retained for every dynamic run. The MAT
snapshot is the M2 import boundary; the instrument exports remain the primary provenance
source described in `docs/data_dictionary.md`.

### Loading

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
