# Recorded-time intervals for advanced Vyazovkin

The audit checks the full measured conversion branch over `[α − Δα, α]` on 18 pinned
pICNIK TG files. It uses the same 113 run and mass-reference scenarios and 15 targets as
the [TG conversion audit](public_tg_conversion.md). The tested increment widths are
0.01, 0.02, and 0.05. The lower target must be positive, so a target equal to the
increment width is outside the admissible interior.

Both boundaries must be unique upward crossings within the selected reaction record.
The method clips the recorded clock and temperature at the two boundaries with
linear interpolation and retains every acquired point between them. An exact
boundary sample appears once. The integration path retains Kelvin temperature in
acquisition order and time in minutes.

Strict branch support requires conversion to rise at every step of the clipped
path. A reversal can make interior levels ambiguous even when the two boundary
crossings are unique. Plateaus likewise have several times at a single conversion
level. The audit retains their source brackets, reversal and plateau counts, largest
conversion backtrack, and total plateau duration. It uses the raw normalization
values, with no noise tolerance or conversion repair. Recorded temperature can
rise or fall along a valid conversion branch. The script reports temperature
reversals, time steps, and measured node counts separately.

## Common recorded-time support

Of 5,085 requested run, scenario, target, and width combinations, 4,388 have strict
support. There are 565 requests outside the conversion interior, 93 with a missing
or ambiguous boundary, 35 with a reversal inside the branch, and 4 with an internal
plateau. The 39 interior reversal and plateau cases end at α = 0.025 or 0.05.
For example, the 1 K/min PP trace normalized to its endpoints has a reversal of
0.0050 in the window ending at α = 0.025 with Δα = 0.02.

The entries below count common targets from the same 15-point grid, in the order
Δα = 0.01, 0.02, and 0.05.

| Material group | Plateau-median references | Intersection across reference and window scenarios |
|---|---:|---:|
| Constant energy | 14 / 14 / 12 | 14 / 14 / 12 |
| Two steps | 14 / 14 / 12 | 11 / 11 / 9 |
| Variable energy | 14 / 14 / 12 | 9 / 9 / 9 |
| Polypropylene | 12 / 12 / 12 | 11 / 11 / 10 |

Every baseline group supports all nine central targets (α = 0.10–0.90) at each
width. That gives 486 supported run-level intervals in the central plateau-median
case. For PP the common set remains the same after omitting the 10 K/min trace.
Across all PP scenarios the supported targets are the nine central points plus
0.95 and 0.975 at widths 0.01 and 0.02. Width 0.05 additionally loses α = 0.10
because the lower target at 0.05 falls outside a cropped reaction record.

All strictly supported intervals span more than one acquisition bracket. The
variable-energy 2.5 and 5 K/min traces nevertheless contribute just one actual
interior point at α = 0.30 and Δα = 0.01 under plateau-median references. The
clipped duration is about 0.310 and 0.165 min, respectively.

## Fixed-energy quadrature comparison

At each supported interval the audit evaluates the time integral of
`exp(−E_probe/(R T(t)))` between the two crossings. It compares trapezoids on the
clipped measured nodes, 41 uniformly spaced conversion nodes as in the legacy
sampling rule, and two time refinements. Each refinement subdivides every clipped
segment into 8 or 16 equal time steps with linear interpolation of its measured
temperatures. The legacy comparison uses the same measured branch so the
difference isolates sampling. All calculations use fixed 60, 120, and 240 kJ/mol
probes.

The report records each log integral and the relative change from the 16-subdivision
value. Agreement under refinement measures numerical error on the piecewise-linear
observed history. Acquisition density and uncertainty in the underlying
temperature history are separate. A supported window contained within a single
acquisition bracket, for example, has only two source temperatures regardless of
the subdivision count.

Across all strictly supported intervals and probes, the largest absolute relative
differences from the 16-subdivision result are 0.0820% for measured-node
trapezoids, 0.398% for legacy conversion resampling, and 0.000965% from 8 to 16
subdivisions. The largest legacy discrepancy occurs in the two-step 10 K/min
trace under the earlier-terminal scenario at α = 0.80, Δα = 0.05, and a 240 kJ/mol
probe. These comparisons address quadrature on the recorded program. The
physical interpretation of each conversion window and the effect of Δα on an
estimated profile require separate analysis.

## Reproduction

```sh
julia --startup-file=no --project=. scripts/audit_public_tg_vyazovkin.jl
```

The [audit settings](../config/public_tg_vyazovkin.toml) and
[machine-readable report](audits/public_tg_vyazovkin.toml) retain the probe settings,
source fingerprints, and per-window branch and quadrature evidence. Each material
group reports the common supported targets for every scenario and width. The
polypropylene group also reports those targets after omitting the run near
10 K/min with negative apparent terminal mass. The Julia audit reproduced byte
for byte, and an independent Python reconstruction matched all 5,085 classifications
and the retained path counts, source brackets, times, and temperatures.
