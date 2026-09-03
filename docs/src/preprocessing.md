# Preprocessing

M3 transforms raw calibration experiments into traceable `ProcessedExperiment` values. The
transformation never mutates an `Experiment`, does not clip conversion, and does not force a
monotonic curve.

## Default workflow

```julia
using IsoconversionalAnalysis

config = load_config("config/analysis_defaults.toml")
catalog = load_dataset_catalog("config/datasets.toml")
experiments = load_experiments(catalog; roles=[:calibration])
processed = preprocess(first(experiments), config)
```

The accepted `m3_default_v1` profile performs these operations in order:

1. Drop non-finite core rows under the explicit `drop_nonfinite` policy.
2. Select dynamic segment 2 when recorded; otherwise retain all finite rows and emit a
   fallback warning.
3. Estimate mass at exactly 150 and 700 °C with robust local-linear fits in ±2 K windows.
4. Retain the inclusive 150–700 °C interval and optionally rebase its time coordinate.
5. Calculate raw conversion without clipping:

   ```math
   \alpha(T) = \frac{m_{150} - m(T)}{m_{150} - m_{700}}.
   ```

6. Fit cubic local polynomials against time using observations in a ±5 K physical
   temperature neighborhood.
7. Calculate local `dalpha_dt_min_inv` and `heating_rate_K_per_min`, then form
   `dalpha_dT_K_inv = (dalpha/dt)/(dT/dt)`.
8. Integrate `dalpha_dT_K_inv` along the recorded temperature path and compare it with the
   configured analysis conversion.

Every retained value carries its original acquisition row through `source_row_indices`.
`config_fingerprint` is a deterministic SHA-256 identifier covering the complete
preprocessing profile and both reference temperatures.

## Diagnostics and failure behavior

`PreprocessingDiagnostics` records row-selection counts, fitted reference masses, measured
heating rate, conversion range, reversal and negative-derivative counts, reconstruction
errors, and warnings. Values slightly or materially outside `[0, 1]` remain visible. The
default monotonic policy diagnoses decreases but never repairs them.

Structural failures throw `PreprocessingError`, including an uncovered reference interval,
non-positive reference mass loss, non-increasing retained time, and a non-positive fitted
heating-rate curve. Setting `invalid_row_policy = "error"` or
`monotonic_conversion_policy = "error"` promotes those configured conditions to failures.

## Sensitivity evidence

The real-data audit compares raw finite differences, finite differences after smoothing, and
local-polynomial half-windows of 2.5, 5, and 10 K on the four 50 wt% calibration runs. It also
compares robust and linear-interpolated mass references and repeats conversion with 120 °C,
650 °C, 750 °C, and the highest common observed temperature.

The tracked evidence files are `docs/preprocessing_audit.md`,
`docs/audits/m3_preprocessing_sensitivity.toml`,
`docs/figures/m3/preprocessing_qc.svg`, and
`docs/figures/m3/preprocessing_sensitivity.svg`.

Regenerate all evidence with:

```sh
julia --startup-file=no --project=docs scripts/audit_preprocessing.jl
```

CairoMakie is isolated in the documentation environment. The numerical package returns
arrays and diagnostics and does not import a plotting backend. The 15 ramp-and-hold
experiments remain held out and did not influence the M3 default.

## Numerical utilities

```@docs
preprocess
preprocessing_fingerprint
finite_difference_derivative
local_polynomial_estimate
cumulative_trapezoid
```
