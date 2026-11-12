# Legacy comparison

M9 compares the corrected Julia workflow with preserved thesis data and records findings
from the earlier review of the MATLAB implementation. The original source is stored separately.

## Numerical comparison boundary

Exact numerical parity is required at the raw dynamic-data boundary. All 20 named arrays in
`code/2021corridas/2021corridas.mat` are exactly equal, including NaNs, to their configured
`St01`–`St20` aliases in `code/2021mar.mat`.

The audit also establishes that every four-column array in `code/2021marcorrected.mat` is an
exact contiguous copy of the corresponding raw array's first four columns. The arrays begin
near 110 °C and end near 794–799 °C. They contain no stored normalized conversion, smoothed
derivative, fit, or uncertainty result, so this file is a cropped raw reference rather than a
preprocessing regression target.

## Corrected-method policy

Beyond raw import, comparison is method-level. The Julia implementation preserves the
scientific questions and named method families while correcting identified defects and
adding explicit units, provenance, diagnostics, sensitivity analysis, synthetic recovery,
and held-out prediction gates. In particular, it does not reproduce the legacy Friedman
variable errors, the hard-coded advanced-Vyazovkin objective target, order-dependent peak
fits, duplicated kinetic intercepts, or workspace-dependent simulations.

Saved result and publication workspaces are retained as qualitative evidence. The archive
does not record an unambiguous producing command, environment, script revision, and variable
lineage for those files, so their activation energies and fitted parameters are not release
acceptance oracles.

The complete comparison table, 13 intentional deviations, 20 per-experiment raw/crop checks,
and evidence-tier policy are in `docs/legacy_comparison_audit.md`. Machine-readable evidence
and SHA-256 hashes for the 18 MAT artifacts used by the comparison are in
`docs/audits/m9_legacy_comparison.toml`.

Regenerate both artifacts from the repository root with:

```sh
julia --startup-file=no --project=. scripts/audit_legacy_comparison.jl
```

Any future loss of exact raw-alias equality is a regression. Differences in preprocessing,
activation energies, peak parameters, reaction-model coefficients, or predictions must be
judged against the corrected Julia equations and their M3–M8 validation gates, not numerical
similarity to an ambiguous legacy workspace.
