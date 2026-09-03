# Configuration and logging

The default configuration is `config/analysis_defaults.toml`. It records the primary
150–700 °C conversion interval selected in M0, the endpoint sensitivity cases, canonical
boundary units, output location, logging threshold, and accepted M3 preprocessing profile.

```julia
using IsoconversionalAnalysis

config = load_config("config/analysis_defaults.toml")

with_project_logger(config) do
    @info "analysis_started" configuration = config.source_path
end
```

Relative paths in a configuration file are resolved relative to that file. The same command
therefore has the same result regardless of the caller's current working directory.

The `[preprocessing]` table makes row handling, segment selection, reference-mass estimation,
smoothing, differentiation, heating-rate fitting, monotonicity policy, downstream conversion
range, and reconstruction tolerance explicit. Temperature windows are half-widths in kelvin,
not sample counts.

The `[isoconversional]` table freezes the common conversion grid, interpolation policy,
confidence level, minimum experiment count, Friedman rate floor, and the R² warning threshold.
Its nested `[isoconversional.advanced_vyazovkin]` table records the conversion interval,
integration resolution, energy bounds, optimizer controls, and boundary/flatness tolerances.

The `[deconvolution]` table records the candidate peak counts and information criterion;
conversion interval and deterministic point budget; height, skew, width, center-separation,
and minimum-area constraints; independent and joint multistart counts; optimizer tolerances;
and the local identifiability, boundary, and confidence thresholds. These settings are part of
each result fingerprint.

The `[reaction_models]` table freezes the ordered candidate registry, the fixed-model subset
used by the compensation diagnostic, exact conversion fit interval, deterministic point and
rate thresholds, BIC/AICc choice, activation-energy and empirical-exponent plausibility
bounds, confidence level, identifiability and ambiguity thresholds, and the
leave-one-heating-rate-out predictive error gate. These settings are also part of each M6
result fingerprint.

Unknown enum values, non-positive windows or tolerances, invalid polynomial degrees, malformed
ranges, inconsistent advanced-Vyazovkin grids, and invalid deconvolution constraints fail
configuration loading with `ConfigurationError`. The same applies to duplicate or unsupported
reaction models, an empirical model in the fixed compensation subset, and invalid kinetic
bounds or thresholds.

Julia's logging macros attach key-value pairs, such as `configuration` above, as structured
metadata. The M1 console backend is intentionally simple and can later be replaced without
changing scientific routines.
