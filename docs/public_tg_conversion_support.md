# Common conversion crossings on public TG data

The audit compares 13 simulated pICNIK traces in the constant-energy, two-step, and
variable-energy groups, and five experimental polypropylene (PP) traces. Each group
uses the same material label and the reaction intervals and mass references from the
[TG conversion audit](public_tg_conversion.md). The 15 tested conversion values run
from 0.01 to 0.99. The central grid contains nine points from 0.10 to 0.90.

A target has common support when every run in the group has exactly one upward
crossing in acquisition order. A linearly interpolated crossing uses two observations
from the recorded segment. The audit retains source line numbers and recorded time and
temperature at each crossing. Missing crossings, repeated crossings, reversals,
plateaus, tangencies, and exact endpoints remain visible as separate cases. The
1,695 per-run evaluations returned 1,649 unique upward crossings, 42 missing targets,
and four ambiguous targets. The public runs had no exact target-valued sample on this
grid. Synthetic tests exercise exact-sample and plateau handling.

## Coverage with plateau median references

| Group | Runs | Common targets | Temperature range at α = 0.50 (K) | Central spans across runs (K) |
|---|---:|---|---:|---:|
| Constant energy | 5 | 15 of 15 | 600.75–684.26 | 67.52–96.73 |
| Two steps | 4 | 15 of 15 | 622.39–712.04 | 70.36–128.56 |
| Variable energy | 4 | 15 of 15 | 474.21–521.47 | 39.12–71.94 |
| Polypropylene | 5 | 14 of 15 | 629.67–708.59 | 65.81–84.40 |

The last column gives the range of temperature spans across the nine central targets.
At α = 0.50 the five PP runs have a reciprocal temperature span of 176.89 × 10⁻⁶ K⁻¹.
The nominal 2.3 K/min PP trace crosses α = 0.01
upward, downward, and upward, so that target lacks common support under the plateau
median references. The [figure](figures/public_tg_conversion_support.svg) shows the
measured temperature coverage at the supported grid points.

## Sensitivity to reaction windows and run selection

The intersection across all reference and reaction-window scenarios retains the
following tested points:

| Group | Targets common to all tested scenarios |
|---|---|
| Constant energy | All 15 |
| Two steps | 12 from α = 0.01 through 0.90 |
| Variable energy | Nine from α = 0.40 through 0.99 |
| Polypropylene | Eleven from α = 0.10 through 0.90, plus 0.95 and 0.975 |

With fixed plateau medians, the 150–800 °C crop removes late targets in the two-step
group and early targets in the variable-energy group. The 300–500 °C crop removes early
PP targets. Imposing a zero terminal reference removes α = 0.99 for PP. This intersection
describes the scenarios tested on the grid. For the two-step group, the earlier
terminal reference retains all 15 targets while shifting a crossing by as much as
232 K. A shared conversion target still depends on the mass-reference model.

Omitting the PP trace near 10 K/min with the negative apparent terminal mass leaves
the same 14 plateau-median targets and the same minimum and maximum temperatures at
each target. It also leaves the 11-target intersection across the seven scenarios.
The omission reduces the centered sum of squares of reciprocal temperatures by
22.57% at α = 0.50 and by 20.60–23.50% across the nine central targets. The extremes
of the observed temperature distribution therefore tell only part of the story. The
audit records `Sxx = Σ(1/Tᵢ − mean(1/T))²` and each leave-one-out profile for all four
groups. Omitting the PP run near 1 K/min reduces the α = 0.50 temperature span from
78.92 to 58.11 K. Omitting the run near 2.3 K/min restores the α = 0.01 crossing
in the plateau-median case. Both span and Sxx describe the observed temperature design.

## Recorded heating programs

The audit measures the slope and temperature residuals of an affine fit to the
recorded temperature history, both across the full record and between the acquisition
brackets at α = 0.10 and 0.90. It also reports the observed temperature/time secants
over α = 0.10–0.50 and 0.50–0.90. The PP run near 20 K/min has corresponding secants
of 19.63 and 20.87 K/min. Its maximum absolute residual to the central-window fit is
0.81 K, compared with 9.95 K across the full record. The five PP traces have strictly
increasing temperatures in their central windows, while four have local decreases
somewhere in the full record. The simulated programs are affine within numerical
precision. The recorded histories can support method-specific applicability checks
before apparent activation-energy estimation.
[Method checks on the recorded histories](public_tg_method_checks.md)

## Reproduction

Run the audit from the repository root with the pinned raw files available.

```sh
julia --startup-file=no --project=. scripts/audit_public_tg_conversion_support.jl
python scripts/verify_public_tg_conversion_support.py
julia --startup-file=no --project=docs scripts/plot_public_tg_conversion_support.jl
```

The [machine-readable evidence](audits/public_tg_conversion_support.toml) records all
25 group and scenario combinations, 24 omission cases, and the crossing-level
provenance. The independent verifier reconstructed all 1,695 per-run classifications,
375 group grid points, and 360 omission grid points from the raw TG values. The
previous conversion audit reproduces byte for byte with the shared scenario builder.
