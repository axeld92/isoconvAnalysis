# Public TG conversion validation

The audit checked 113 mass-reference and reaction-interval cases on 13 simulated
pICNIK TG traces and five experimental polypropylene (PP) traces. It also checked
original-sample mass loss on two mixchar TG examples. The [protocol](../config/public_tg_conversion.toml)
records the windows chosen after inspecting the traces, and the [evidence](audits/public_tg_conversion.toml)
records the actual retained rows, source lines, mass references, conversion extrema,
and crossing counts.

## Conversion and audit procedure

For measured remaining mass \(m\), `tg_conversion` calculates
\(\alpha=(m_0-m)/(m_0-m_f)\) using explicitly supplied initial and terminal references.
Each reference can be a fixed mass in the signal's native unit or the median of a
segment-local row interval. The reaction interval is supplied independently in acquisition
order. The public audit translates the configured temperature windows to rows only
when the matched rows form a contiguous interval.

| Case | Reaction interval | Initial and terminal mass references |
|---|---|---|
| Record endpoints | Full record | First and last observations |
| Plateau medians | Full record | Early and late window medians used for comparison |
| Reaction window | 150–800 °C simulated or 300–500 °C PP | Same as plateau medians |
| Later initial | Full record | Later initial window and original terminal window |
| Earlier terminal | Full record | Original initial window and earlier terminal window |
| Window and references | Same as reaction window | Both shifted windows |
| Zero terminal (PP) | Full record | Early window median and assumed zero terminal mass |

The implementation retains the recorded temperature, time, mass, and source lines.
The audit checks mass reconstruction and equality on overlapping rows when only the
interval changes. It reports progress excursions and reversals, and counts upward
crossings on a 15-point grid from 0.01 to 0.99. An upward bracket has a left value below
the target and a right value at or above it. The first interval also counts when it starts
exactly at the target and rises. Temperature shifts use acquisition-order linear
interpolation where both cases have one bracket at the same conversion. Missing and
multiple brackets remain in the results. Independent checking found zero multirow
plateaus equal to any tested target. A general isoconversional sampler will need to
handle such plateaus explicitly. These window choices and sensitivity ranges are
conditional on the downloaded traces.

## Observed sensitivity

- **PP.** In the 300–500 °C window, progress already ranges from 0.022 to 0.076 at
  the first retained observation. The full comparison window supports a unique
  crossing in every run at each tested point from 0.025 to 0.99. The 300–500 °C
  interval supports that only from 0.10 to 0.99. First/last masses versus plateau
  medians shift central crossings (0.10–0.90) by at most 0.17 K. A later initial
  reference (130–170 °C) shifts central crossings by up to 3.27 K and the 0.01
  crossing by up to 103 K. Low-conversion crossings are sensitive to small early
  mass changes ahead of the main loss.
- **PP terminal mass.** The 10 K/min trace has an apparent median terminal mass of
  −0.1241 mg at 570–590 °C. The other four range from 0.0867 to 0.2159 mg. Using
  zero terminal mass as a sensitivity case shifts central crossings by up to 1.38 K
  and crossings across the full grid by up to 5.75 K. Its last progress readings
  span 0.981–1.012 across PP, including above-one values in the 10 K/min trace.
  These readings require balance-baseline or blank provenance for physical residue
  estimates.
- **Simulated traces.** The 150–800 °C interval retains all 15 grid points at all five
  constant-energy rates. For the variable-energy case, the same interval starts
  as late as \(\alpha=0.349\) and only supports unique crossings across all rates from
  0.40 onward on this grid. The fastest two-step trace reaches just 0.932 by 800 °C.
  Using the 710–750 °C window for the terminal mass while retaining the full record
  produces last values up to 1.214 in the two-step traces and shifts a crossing by up
  to 232 K. Moving the initial reference to 130–170 °C shifts a crossing by up to
  152 K in the variable-energy traces, where some loss is already underway.
- **mixchar.** The reported sample masses of 18.96 and 15.29 mg give terminal loss
  fractions of 0.7630 and 0.7191 of the original sample mass for Juncus and Marsilea.
  The late 775–795 °C medians give 0.7621 and 0.7180. Each file has one temperature
  trace and no recorded time. These are observed whole-sample loss fractions.
  Reaction completion and time-dependent conversion remain undetermined.
  [Pinned mixchar vignette](https://github.com/smwindecker/mixchar/blob/8e26d64d6603e583e6cbcddb1ea6123cd130963f/vignettes/mixchar.Rmd)

The pICNIK publication identifies PP as an experimental case. A companion symposium
abstract reports four PP rates (1, 2.3, 5, and 20 K/min), an SDT Q600 and nitrogen at
100 mL/min. File-level linkage, the additional 10 K/min trace, blank correction, and
replicates remain undocumented in that abstract. The raw files supply temperature,
time, and reported mass. [Publication record](https://doi.org/10.1016/j.cpc.2022.108416)
· [Symposium proceedings, pp. 105–106](https://politecnicojic.edu.co/images/downloads/investigacion/grupos/camer/productos/iii-simposio-internacional-camer-2022.pdf#page=111)

## Reproduction

Download the pinned pICNIK and mixchar inputs. The full conversion audit also
uses the optional `rdata` decoder to prepare the two mixchar CSVs.

```sh
python -B scripts/fetch_public_benchmarks.py --sources picnik mixchar
python -m pip install --only-binary=:all: rdata==1.0.0
python -B scripts/prepare_mixchar_examples.py
julia --startup-file=no --project=. scripts/audit_public_tg_conversion.jl
```

The command regenerates the TOML evidence with pinned source checksums. The package
tests cover fixed and observed references, acquisition rows, segment selection,
absent time, invalid references, mass scaling, and apparent negative residue. An
independent Python CSV calculation agreed with all 113 cases and 1,695 crossing
evaluations.
