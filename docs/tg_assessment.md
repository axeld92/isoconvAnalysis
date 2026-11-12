# Configurable ICTAC-aligned TG assessment

Run the assessment from the repository root with the Julia project instantiated.
Pass a dataset config and an output directory. Each directory receives
`tg_assessment.md` and `tg_assessment.toml`.

Julia 1.13 uses the version-specific lockfiles. They include the SHA 1 standard
library and updated RuntimeGeneratedFunctions and PrecompileTools for
precompilation on this Julia series.

On Julia 1.13.3, both datasets reproduced the archived activation-energy profiles
exactly. Replaying the archived frozen pine fit also reproduced the hold RMSEs
exactly. Recalibration on 1.13.3 changed the hold RMSEs by at most 0.0033 percentage
points. The cross-version refit regression allows 0.0001 mg, one increment of the
source balance. The archived-fit replay uses a separate, tighter numerical check.

```sh
julia --startup-file=no --project=. scripts/assess_tg.jl config/tg_assessment_pine.toml results/tg_assessment/pine
julia --startup-file=no --project=. scripts/assess_tg.jl config/tg_assessment_polypropylene.toml results/tg_assessment/polypropylene
```

For pine, use the extracted TG files and acquisition manifest described in the
[white-pine benchmark](public_pine_tg.md). For polypropylene, fetch the pinned
inputs with `julia --startup-file=no --project=. scripts/fetch_public_tg.jl`.
The two example configs reproduce the existing conversion reference choices
and method settings. The [polypropylene benchmark](public_tg_energy.md)
describes the source and known-truth controls.

After generating both reports, run the focused regression in PowerShell.
It compares the new reports with the archived benchmark results.

```powershell
$env:ISOCONV_TG_ASSESSMENT_RESULTS = 'results/tg_assessment'
'using Test; include("test/test_tg_assessment.jl")' | julia --startup-file=no --project=. -
```

## Configure a dataset

Use the [pine config](../config/tg_assessment_pine.toml) or the
[polypropylene config](../config/tg_assessment_polypropylene.toml) as a starting point.

- Set `schema_version = 1`, `id` and `label`. Declare the time, temperature, and
  remaining-mass columns and their units in `[columns]`. The default delimiter
  is a comma. Set `delimiter = "\t"` for tab-separated exports and set `segment`
  to the segment column name when one is present. The importer also accepts
  `encoding` and `header_prefix` for instrument exports that require them.
- Declare inputs in `[[runs]]` with `id` and `path`, or point `manifest` to an
  existing TOML inventory. `manifest_table` defaults to `runs`. An input can
  also declare `sha256`, `rows`, `selected_segment` and `sample_mass_mg`.
  `manifest` and `data_root` resolve relative to the config directory. Run
  paths resolve under `data_root`, which defaults to the config directory.
  For a file with multiple segment labels, `selected_segment` identifies one
  segment by its original ordinal. The selected segment supplies the time and
  temperature history for that run.
- Put all repeats of each thermal program in one `[[programs]]` entry. List
  their IDs under `runs` and set `kind` to `ramp` or `hold`. Assign one `role`
  from `calibration`, `mass_series` and `holdout`. Set `expected_repeats` if
  the group size is known. Specify `nominal_rate_K_per_min` when known and
  `nominal_mass_mg` for each sample-mass group. The assessment can use measured
  rates for ramp programs without a declared nominal rate. When withholding a
  ramp, declare nominal rates for that ramp and every calibration ramp so the
  command can check for overlap. Holds use the `holdout` role. Use
  `measured_ramp_window_K` to limit measured rate estimates to the ramp.
- In `[conversion]`, declare `initial` and `terminal` mass references and
  optionally a `reaction` interval. Select recorded windows with
  `temperature_K` or `temperature_C`, `rows`, or `last_seconds`. A reference
  can also be a declared `mass` in the source unit. A reaction interval can
  start at the first observed sample at or above `start_K`. Use
  `[conversion.variants.<name>]` to compare alternate choices. Set
  `repeatability_alpha` to select the conversion targets for repeat summaries.
  The optional `omit_run` adds a leave-one-run-out energy sensitivity case.
- Under `[evidence]`, document the experimental source and give each available
  check a `status` and `note`. Use `metadata`, `calibration`, `repeatability`,
  `prediction` and `components` for the experimental checks. The command
  derives the repeatability status from the roster and the executed prediction
  status from the scores. The reports show missing evidence and preserve
  individual method statuses at unsupported conversion targets.

## Predict an excluded program

When a `[prediction]` block is present, the command fits a constant-energy
mass-loss baseline once using the calibration ramps. It freezes the fit before
scoring every run in the holdout programs. The mass-series role contributes
sample-mass diagnostics. For the pine config, the fit uses 18 main ramps and
the 593 K program contributes all three holdout traces.

`start_K` can set the conditioning point for each measured program. Each
prediction uses that first mass and the recorded time and temperature history.
`training_window` declares the additional fit window. `holdout_window` selects
the evaluation interval. The pine config uses the final contiguous interval
within 590–596 K. Optional `checkpoints_min` and `final_window_min` provide
additional checks on the hold. Prediction requires a mass unit convertible to
mg. Mass-percent traces can be assessed using conversion and energy diagnostics.
Use `settings_file` and `[prediction.settings]` to supply baseline fit settings
and overrides as in the pine config.

The TOML records input hashes, conversion references, method and run statuses,
repeat summaries, and the calibration and holdout IDs. For a prediction it
also records the fitted parameters and time-weighted errors on the evaluation
window and the full record after conditioning. An ICTAC status describes the
evidence for an individual recommendation. Review the calibration, measurement
and uncertainty evidence alongside the numerical results.
