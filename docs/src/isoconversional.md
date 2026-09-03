# Isoconversional analysis

M4 estimates apparent activation-energy profiles from groups of processed dynamic
experiments at a common conversion. The package implements Friedman, KAS, FWO, Starink,
and advanced Vyazovkin behind one typed API.

```julia
using IsoconversionalAnalysis

config = load_config("config/analysis_defaults.toml")
catalog = load_dataset_catalog("config/datasets.toml")
raw = load_experiments(catalog; roles=[:calibration])
processed = preprocess.(raw, Ref(config))

group = filter(experiment -> experiment.composition.waste_tire == 0.5, processed)
starink = analyze(:starink, group, config)
advanced = analyze(:advanced_vyazovkin, group, config)
```

`IsoconversionalResult` retains every requested conversion. A point that cannot be evaluated
has a missing estimate and an explicit `:invalid` status; warnings and optimizer-boundary
solutions remain visible in `point_diagnostics`. Input experiments are canonically ordered
by identifier, making numerical output and provenance fingerprints independent of caller
order.

## Shared interpolation contract

All methods use acquisition-order, first-upward-crossing linear interpolation. The policy
does not extrapolate, clip conversion, or silently make a measured curve monotonic. If noise
creates multiple upward crossings, the earliest is used and the crossing count is reported.

## Regression methods

Friedman regresses `log(dalpha/dt)` against `1/T`. KAS, FWO, and Starink use their documented
integral approximations with the measured global heating rate. Their confidence limits are
two-sided Student-t intervals for the ordinary-least-squares slope, transformed to
activation energy. R² is a fit diagnostic, not a confidence level; points below the
configured warning threshold remain available with `:warning` status.

FWO is retained for comparison with the thesis and literature. Starink is preferred when a
linear integral approximation is wanted because it has lower approximation error over the
usual kinetic range.

## Advanced Vyazovkin

The advanced method integrates each experiment over a conversion interval, evaluates the
full pairwise objective for any experiment count, and minimizes that objective directly.
Arrhenius integrals are evaluated in a log-scaled form. A bounded scan brackets a
golden-section refinement and records boundary, flatness, and multiple-minimum diagnostics.

The reported confidence interval follows the Vyazovkin–Wight pairwise-variance Fisher
procedure. It is intentionally distinct from the objective minimum and from conversion-
endpoint sensitivity.

## Validation evidence

The synthetic fixtures enforce known-energy recovery, noisy-data tolerances, arbitrary
experiment counts, confidence-root evaluation, visible boundary solutions, and permutation
invariance. The real-data evidence in `docs/isoconversional_audit.md` reports all five methods
for five compositions, plus conversion-endpoint sensitivity. The 15 ramp-and-hold runs are
not used in M4.
