# Kinetic triplets and reaction models

M6 fits an overall-conversion Arrhenius triplet jointly to the four dynamic heating programs
for one composition. It deliberately does not propagate the M5 peaks: their cross-rate
identities were not established.

## Model registry and units

The public registry contains conventional differential equations rather than mechanistic
assignments:

```julia
models = reaction_model_registry()
value = reaction_function(:sestak_berggren_2, 0.4; m=0.25, n=1.1)
```

The fixed candidates are reaction orders `f1`–`f4`, Avrami–Erofeev `a2`–`a4`, contracting
geometry `r2`–`r3`, diffusion `d1`–`d4`, power laws `p2`–`p4`, and `random_scission`.
The two empirical candidates are

```text
sestak_berggren_2: f(alpha) = alpha^m * (1-alpha)^n
sestak_berggren_3: f(alpha) = alpha^m * (1-alpha)^n * [-ln(1-alpha)]^p
```

All functions require `0 < alpha < 1` and return a positive dimensionless value. Public
activation energies use kJ/mol. Because time is measured in minutes, `A` uses `min^-1`.

## Joint multi-rate fit

For a fixed reaction function, M6 fits

```text
ln(dalpha/dt) - ln[f(alpha)] = ln(A/min^-1) - E/(R*T)
```

with one intercept only. The empirical models add `ln(alpha)`, `ln(1-alpha)`, and optionally
`ln[-ln(1-alpha)]` as linear design columns. This removes the duplicated intercept and
pre-exponential terms present in the legacy scripts.

```julia
config = load_config("config/analysis_defaults.toml")
fit = fit_kinetic_triplet(composition_group, :sestak_berggren_2, config)
comparison = compare_reaction_models(composition_group, config)
selected = only(filter(result -> result.model.name == comparison.selected_model,
                       comparison.results))
```

The fit interval uses exact first-upward conversion crossings. Each experiment contributes a
deterministically bounded number of points. Nonpositive observed rates are counted and
excluded because logarithms are undefined; they are never replaced with an arbitrary floor.

`compare_reaction_models` retains all candidates, the raw BIC/AICc minimum, structural
eligibility, and the selected model. Full design rank, configured condition number, activation
energy, and empirical-exponent bounds determine structural eligibility. A criterion tie and
failure of the leave-one-heating-rate-out error gate remain explicit statuses.

## Uncertainty and prediction

`KineticTripletResult` includes coefficient estimates, covariance and correlation matrices,
local two-sided Student-t intervals, residual metrics, Durbin–Watson statistics, information
criteria, and every held-out-heating-rate result. These intervals are conditional on the
chosen model and treat serially correlated derivative samples as independent, so the
heating-rate-level validation is the stronger predictive check.

```julia
rate = predict_rate(selected, 700.0, 0.5)
lower, central, upper = predict_rate_confidence_interval(selected, 700.0, 0.5)
```

The confidence interval is for the fitted mean log-rate transformed back to rate space. It is
not a prediction interval for a future measurement and does not include model-selection or
preprocessing uncertainty.

## Compensation-effect diagnostic

```julia
diagnostic = analyze_compensation(first(composition_group), 180.0, config)
```

This fits the fixed registry to one heating program, regresses the resulting `ln(A)` values
on `E`, and evaluates the line at an independently supplied reference energy. The result
includes a mean-response interval. The parameter pairs share the same data and their apparent
linear correlation can arise mathematically, so a high R² is not evidence of a physical
compensation law and the estimate is not substituted for the direct multi-rate fit.

## M6 real-data decision

The complete audit is in `docs/reaction_model_audit.md`, with all 95 candidate fits, 380
heating-rate validation rows, and 20 compensation diagnostics in
`docs/audits/m6_reaction_model_results.toml`. Every selected constant-triplet model fails the
predeclared held-out-rate log-RMSE gate. M7 will therefore retain these fits only as descriptive
baselines and evaluate a conversion-dependent isoconversional description for prediction.

