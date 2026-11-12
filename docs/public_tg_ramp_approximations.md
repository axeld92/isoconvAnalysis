# Constant-ramp checks on public TG data

KAS, FWO and Starink replace the time integral for a linear heating program with a
temperature integral divided by the constant heating rate. They then use distinct
approximations of the temperature integral. The ramp audit checks these two
assumptions separately at fixed probe energies of 60, 120 and 240 kJ/mol.
[Starink's comparison of linear-heating methods](https://eprints.soton.ac.uk/18822/)

The audit uses all 18 pinned pICNIK TG runs and the reaction intervals and references
from the [public TG conversion audit](public_tg_conversion.md). A prefix starts at the
first retained row and ends at the unique upward crossing. The plateau-median case
uses all 15 conversion targets and also tests a later start at the first measured
temperature at or above 150 °C for synthetic TG or 300 °C for polypropylene (PP).
The other reference and crop scenarios use the nine central targets from 0.10 through
0.90. The report carries both the full PP group and an omission of the nominal
10 K/min run.

## Recorded program and finite start

The comparison ramp shares the observed start temperature, crossing temperature
and elapsed time. The audit reports the endpoint secant, the sample-weighted OLS
slope, two half-duration secants, the range of step slopes, cooling and flat steps,
and temperature residuals from the comparison ramp. The numerical integral
`J = ∫ exp(−Eprobe/(R T(t))) dt` follows the piecewise-linear recorded temperatures
and the measured clock. A second integral replaces only the temperature path
with the endpoint-matched affine ramp. Their log ratio isolates the path effect
under these matched endpoints. The slope relative to the full-record OLS slope
is reported separately.

| Material group | Plateau-median common central targets | Maximum absolute prefix temperature residual | Maximum absolute log integral ratio |
|---|---:|---:|---:|
| Constant-energy synthetic | 9 of 9 | < 3 × 10⁻¹³ K | < 2 × 10⁻¹⁴ |
| Two-step synthetic | 9 of 9 | < 3 × 10⁻¹³ K | < 2 × 10⁻¹⁴ |
| Variable-energy synthetic | 9 of 9 | < 3 × 10⁻¹³ K | < 2 × 10⁻¹⁴ |
| Polypropylene | 9 of 9 | 11.78 K | 0.04569 |

The maxima use the first recorded start and all three probe energies. The three
synthetic groups give linear-ramp controls. The largest PP discrepancy comes from
the 20 K/min run at α = 0.90 and 240 kJ/mol. Its measured integral is 4.47% below
the endpoint-matched ramp integral. Omitting
the 10 K/min PP run leaves this worst case unchanged. With all PP reference and
crop scenarios, the maximum absolute ratio at central targets is 0.04625 on the
log scale.

Starting PP at the first measured row at or above 300 °C reduces the maximum
temperature residual to 1.06 K and the maximum absolute log integral ratio to
0.02647. This later start also omits up to **80.5%, 66.6% and 45.6%** of the
measured integral leading to α = 0.10 for the 60, 120 and 240 kJ/mol probes,
respectively. Omitting that history changes the integral even when the retained
temperature program is close to affine. For each evaluated prefix, the evidence
also reports the fraction of the ideal 0 K-to-crossing ramp integral below the
selected start. That fraction and the discarded measured integral are separate
quantities. [Starink's finite-start analysis](https://eprints.soton.ac.uk/19255/)

## Approximation shape

For `x = Eprobe/(RT)`, the numerical reference is
`p(x) = ∫ₓ^∞ exp(−u)/u² du`. The exact affine time integral is
`(Eprobe/R)/β × [p(x_crossing) − p(x_start)]`. The approximation comparison uses
the zero-temperature integral shapes:

| Method | Approximate log shape up to a constant |
|---|---|
| KAS | `−x − 2 log(x)` |
| FWO | `−1.052 x` |
| Starink | `−1.0008 x − 1.92 log(x)` |

The natural-log coefficients match the package's linear-method implementation.
The comparison uses the spread across runs of the approximate log shape minus
`log(p(x))` at the same conversion and probe energy. An additive constant in an
approximation cancels in this spread. The FWO evidence flags whether all crossing
`x` values lie in the interval 20–60 discussed by Flynn and Wall for Doyle's
approximation. [Flynn and Wall's original treatment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6624709/)
· [Starink's linear-heating comparison](https://eprints.soton.ac.uk/18822/)

At α = 0.50 and the 120 kJ/mol probe, the PP shape spreads are **0.00931** for
KAS, **0.09411** for FWO and **0.00190** for Starink. The recorded-versus-affine
log integral ratios span 0.01433 at this point. These values describe separate
sources of variation in the integral. Energy accuracy depends on the fitted
model and the range of conditions in the data. PP is suitable for a conditional
constant-rate comparison with the observed history and formula effects tracked
in the energy-estimation step. Method admission requires a tolerance tied to
the intended use.

## Reproduction

For a prepared TG conversion curve, evaluate a single prefix by supplying the
starting row within its selected segment.

```julia
diagnostic = tg_ramp_diagnostic(
    curve, 0.50;
    start_row=first(curve.reaction_rows),
    probe_energies_kJ_per_mol=[60.0, 120.0, 240.0],
)
```

```powershell
julia --startup-file=no --project=. scripts/audit_public_tg_ramp_approximations.jl
```

The [machine-readable evidence](audits/public_tg_ramp_approximations.toml) contains
1,395 run-prefix checks, 1,363 evaluated prefixes and 4,089 integral probes.
The 4- and 8-node rules applied to each recorded piece differ by at most
1.43 × 10⁻¹⁴ on the log scale on these data. An independent exponential-integral
calculation checked the affine integrals and formula shapes at all 4,089 probes.
The affine log integrals agree within 2.85 × 10⁻¹⁴. These numerical checks address
quadrature over the interpolated path. Temperature and conversion-reference
uncertainty remain separate.
