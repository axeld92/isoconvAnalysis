# Local conversion rates for Friedman analysis

The [TG conversion support audit](public_tg_conversion_support.md) identifies unique
raw upward crossings. `tg_friedman_rate` evaluates the local conversion derivative at
each of those crossing times on the recorded clock. A Gaussian-weighted polynomial
is fitted to the observations inside a fixed, two-sided window. The rate is the
coefficient of the linear term after accounting for the time scale, in min⁻¹. Each
side needs at least two observations and the total needs at least degree plus two.
The audit retains signed rates, fit residuals, observed conversion ranges and the source
lines used by the fit. Missing crossings, time, or sampling support have separate
statuses.

The half-windows are 1.5%, 3%, and 6% of the elapsed time between the plateau-median
crossings at α = 0.10 and 0.90 for each run. These widths scale with the observed
duration of the event across programs. They stay fixed in minutes when the mass
references or reaction interval change. All three widths use degree two. A degree-one
fit at the 3% half-window isolates sensitivity to polynomial degree. The widths and
degrees define comparison scenarios. A second, slower mass-loss region can lengthen
the event duration and make these windows broader than the local feature being sampled.

The log-rate span is the difference between the largest and smallest log of the
positive rate in min⁻¹ across the four settings. At each conversion, the audit
also computes the maximum span of each pair's log-rate difference across settings.
That second quantity measures how derivative settings change the contrast between runs.
Positivity establishes a defined log response at the stated settings. The span and
fit evidence are descriptive and depend on the selected observations and windows.

## Results

The audit covers 18 pICNIK runs, 113 conversion scenarios at the run level, 15
target conversions and four derivative settings. Of the 6,780 evaluations, 6,592
returned positive rates, 184 lacked a unique upward crossing, and four lacked
two-sided sampling. All four sampling failures came from one crossing in the
variable-energy 2.5 K/min run at α = 0.01 after applying both the reaction crop and
the changed references. Only one recorded observation preceded the crossing within
each window.

| Group | Plateau-median targets positive across all four settings | Largest central within-run log-rate span | Largest central change in pairwise log-rate contrast | Central observations per window |
|---|---:|---:|---:|---:|
| Constant energy | 15 of 15 | 0.00462 | 0.00553 | 6–34 |
| Two steps | 15 of 15 | 0.34734 | 0.14083 | 23–107 |
| Variable energy | 15 of 15 | 0.27940 | 0.02994 | 16–73 |
| Polypropylene | 14 of 15 | 0.02068 | 0.01444 | 15–281 |

The central summaries use the nine tested points from α = 0.10 through 0.90. In the
two-step and variable-energy controls, the widest central windows span up to 0.294
and 0.303 in conversion, respectively. The largest within-run spans represent
rate ratios of about 1.42 at α = 0.70 for the two-step control and 1.32 at α = 0.40
for the variable-energy control. Changing the variable-energy initial mass reference
to `later_initial` raises the largest central pairwise contrast span to 0.816 at
α = 0.30, even though all four settings still return positive rates there.

Polypropylene retains all nine central targets under every mass-reference and
reaction-window scenario. The central plateau-median windows span 0.00562–0.18552
in conversion. The largest central within-run rate ratio is 1.021 for the
plateau-median case and 1.038 for the zero-terminal case. Omitting the 10 K/min trace
preserves the common positive targets and largest central within-run span in each
scenario. At α = 0.99 under plateau medians, all settings give positive rates, and the
largest ratio rises to 1.35. Holding the derivative at reference quadratic, changing
to the later initial or zero-terminal mass reference shifts the central PP log-rate
by as much as 0.135 or 0.102 respectively. Those changes include movement of the
raw crossing at the same reported α under the new references. Selection of a
derivative procedure for kinetic inference remains a separate comparison against
known-truth controls.

## Reproduction

Run from the repository root with the pinned raw pICNIK files available.

```sh
julia --startup-file=no --project=. scripts/audit_public_tg_friedman.jl
```

The [protocol](../config/public_tg_friedman.toml) fixes the settings and the PP omission
case. The [audit TOML](audits/public_tg_friedman.toml) records per-run, per-scenario rate
and provenance evidence and group-level positivity and sensitivity on the original
conversion grid. A second audit execution reproduced the evidence byte for byte.
An independent Python/NumPy computation from the pinned raw files matched all 6,780
classifications and the 6,592 fitted rates to a maximum relative difference of
2.6 × 10⁻¹⁴.
