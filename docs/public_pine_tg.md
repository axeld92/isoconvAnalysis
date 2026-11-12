# White-pine TG software benchmark

The [Braunschweig deposit](https://doi.org/10.24355/dbbs.084-202605201212-0) contains
the nitrogen TG runs contributed to the [MaCFP-4 workshop](https://github.com/MaCFP/macfp-db/releases/tag/macfp-4.0).
The supplied 991,278-byte ZIP passed its CRC checks, and all 46 extracted members
matched their archived bytes. All 33 TG exports differ from the pinned GitHub
snapshot only in line endings. The 134,411 time, temperature and mass rows agree
exactly, so the two copies represent the same TG measurements. The [acquisition
manifest](../data/external/macfp_wood_tubs_v2/acquisition_manifest.toml) records a
local SHA-256 for the archive and each member. A publisher checksum was unavailable.
The 12 cone, gasification and microscale combustion CSVs also match the pinned
Git blobs when accounting for line endings. The README changes contributor and
non-TG documentation. Its TG conditions, sample masses and run assignments agree.

## TG conversion and experimental evidence

The 18 main ramps cover 5, 10, 20, 30, 40 and 50 K/min with three repeats each.
Their measured slopes over 450–850 K span 5.00–50.37 K/min. The generic Julia
importer preserved time in seconds, temperature in kelvin and mass in mg. The
exports have no segment column. Conversion uses the observed trajectory from its
first sample at or above 410 K, the median mass over 400–420 K as the post-drying
reference, and the median over the final 60 seconds as the retained mass reference.
The deposited program describes the initial hold, ramp and terminal hold.
Within each main rate, the sample standard deviation of the temperature at
α = 0.50 was 0.12–0.89 K across the three repeats.

| ICTAC assessment area | Result |
| --- | --- |
| Distinct programs and rate spread | Six distinct programs across a nominal decade. |
| Conversion and apparent activation energy | Mass references and reaction rows recorded. Friedman estimated 18 of 19 grid points and recorded-time advanced Vyazovkin estimated 17. Both flagged the common crossing at α = 0.90. Advanced Vyazovkin also flagged its interval at α = 0.95. |
| Acquisition and calibration | Instrument, sample masses, geometry, atmosphere, crucible, calibration materials and baseline correction documented. Blank records and independent thermal-lag evidence are unavailable. |
| Repeatability and mass sensitivity | Three repeats per main rate and a controlled four-mass series at 50 K/min. Mass independence across the full rate range remains unestablished. |
| Predictive validation | The frozen baseline was scored on all three excluded 593 K traces. |
| Uncertainty and step assignment | Reference-window and derivative sensitivities were evaluated. Whole-pipeline uncertainty and the identities of contributing reaction steps remain unresolved. |

At α = 0.50 the apparent energies were 167.1 kJ/mol by Friedman and 167.7 kJ/mol
by advanced Vyazovkin. The common central points span 166–185 kJ/mol. At α = 0.90,
one 5 K/min repeat crosses more than once. At α = 0.95, Friedman rises to 365 kJ/mol
and two 5 K/min repeats reverse within the advanced Vyazovkin interval. Treat that
edge estimate as a diagnostic result with limited support. On common central
points, changing the initial or terminal
reference shifted Friedman by at most 3.8 kJ/mol and advanced Vyazovkin by at most
1.0 kJ/mol. The two alternate derivative windows shifted Friedman by at most
5.1 kJ/mol. These values describe the sensitivity of the calculation to the
declared choices. The numerical audit retains unsupported points and per-run
diagnostics for checking the tail.

### Controlled sample-mass series

Each group is a separate set of three runs at 50 K/min. Residue is expressed as
a percentage of the post-drying reference. Temperature entries are mean ± sample
standard deviation across the three runs.

| Nominal mass (mg) | T at α = 0.10 (K) | T at α = 0.50 (K) | T at α = 0.90 (K) | Residue (%) |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 | 563.19 ± 1.62 | 654.41 ± 0.26 | 702.38 ± 0.86 | 16.77 ± 0.83 |
| 1 | 575.21 ± 0.12 | 655.14 ± 0.11 | 703.84 ± 2.84 | 17.51 ± 1.02 |
| 2 | 566.91 ± 0.19 | 653.77 ± 0.47 | 704.00 ± 0.47 | 16.25 ± 0.67 |
| 4 | 572.02 ± 0.66 | 654.54 ± 0.40 | 713.31 ± 5.53 | 14.89 ± 0.69 |

The midpoint shifts only 1.37 K across group means. At α = 0.10 the range is
12.01 K and the trend with sample mass is nonmonotonic. The 4 mg group crosses
α = 0.90 about 9.5 K later and retains 2.62 percentage points less mass on
average than the dedicated 1 mg group. Those changes support keeping the
sample-mass diagnostic at several conversions and reporting residue separately.

## Frozen prediction on the 593 K program

The benchmark reuses the measured-program, constant-energy nth-order baseline
from the [polypropylene benchmark script](../scripts/audit_public_pp_prediction.jl). The fit uses all 18
main ramps. The 12 mass-series runs and the three holds are excluded from
calibration. The objective gives equal weight to each run and equal weight to
its full post-conditioning record and its 410–760 K window. For each held-out
trace, the solver receives the measured temperature history and the first mass
at or above 410 K. It predicts mass from that observation to the end of the run.
The score for the hold starts at the beginning of the final contiguous 590–596 K
interval, which lasts about 354.5 minutes.

| Hold repeat | Hold RMSE (% of conditioning mass) | Full post-conditioning RMSE (%) | Terminal mass measured (mg) | Terminal mass predicted (mg) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 2.06 | 2.11 | 0.23452 | 0.22006 |
| 2 | 3.89 | 3.89 | 0.25397 | 0.21402 |
| 3 | 4.14 | 4.14 | 0.25028 | 0.20834 |

The baseline predicts too much mass loss by ten minutes into each hold and
underpredicts the terminal mass by 0.014–0.042 mg. The residuals therefore expose
early and late errors even when the whole-hold RMSE is modest. All six training
starts and the refinement converged. Tighter integration changed any predicted
mass by at most 0.000002 mg. The fitted baseline parameters are conditional on
this model and its training programs.

## Reproduction

Use the extracted files and GitHub counterparts at the paths in the
[manifest](../data/external/macfp_wood_tubs_v2/acquisition_manifest.toml) and run
from the repository root with the Julia project instantiated. The benchmark can
run from the verified extracted files alone. It also rechecks the ZIP hash if the
original archive is present at the manifest path.

```sh
julia --startup-file=no --project=. scripts/audit_public_pine_tg.jl
```

The [protocol](../config/public_pine_tg.toml) defines the groups, reference windows
and prediction settings. The [machine-readable audit](audits/public_pine_tg.toml)
records the input hashes, repeat summaries, ICTAC evidence, per-point method
diagnostics, prediction scores and checkpoints.
