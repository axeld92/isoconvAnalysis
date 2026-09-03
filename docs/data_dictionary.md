# Data Dictionary

Last updated: 2026-08-17  
M2 status: Import boundary implemented; remaining provenance gaps are non-blocking  
Scope: legacy dynamic and isothermal thermogravimetric data

## 1. Evidence policy

This document separates three levels of information:

- **Confirmed** — stated by an instrument-export header or demonstrated by exact array comparison.
- **Corroborated** — consistently implied by filenames, values, and multiple legacy scripts, but not present in the instrument metadata.
- **Unresolved** — requires the thesis, laboratory records, or author confirmation.

The Julia implementation must not promote corroborated or unresolved information to fact without recording a decision.

## 2. Canonical dynamic dataset

### 2.1 Source hierarchy

1. **Primary raw source:** the 20 NETZSCH text exports in `code/2021corridas/ExpDat_*wt*hr*.txt`.
2. **Canonical MATLAB snapshot:** `code/2021corridas/2021corridas.mat`.
3. **Legacy alias snapshot:** `code/2021mar.mat`, in which the same arrays were renamed `St01`–`St20`.
4. **Derived comparison artifacts:** `code/2021marcorrected.mat`, `code/2021apr.mat`, `code/corridasjulio.mat`, `code/julio2.mat`, `code/a*neum.mat`, and `code/*_nuevo.mat`.

The text exports are authoritative because they retain instrument metadata and column headings. The canonical MAT snapshot is useful for regression and MAT-import tests, but its arrays contain one trailing all-`NaN` row added during MATLAB import.

SHA-256 fingerprints:

| File | SHA-256 |
|---|---|
| `code/2021corridas/2021corridas.mat` | `78180b021e8b220fb1b8fede20a41467701f3dfd759991c01ecbaed92aac0c82` |
| `code/2021mar.mat` | `aea7fdcd8bde7289a3a552de1dc07235ef2befab27b9e82bd96a46eb23dc619e` |
| `code/isotherm.mat` | `08c0eb9c82691514743061fb44dbbd477e233c36b61e2febddd6a16170032ca8` |
| `code/corridasjulio.mat` | `8b666bb8c6d8fc218a8bf8f9ffc3bfa7102fa3f8cc1ad84cebf6d1e8ed4db821` |

### 2.2 Experimental design

The canonical dynamic set contains 20 distinct experimental conditions:

- nominal waste-tire contents: 0, 25, 50, 75, and 100 wt%
- nominal heating rates: 5, 10, 15, and 20 K/min
- one retained export per composition/heating-rate condition
- nominal temperature program: approximately 25 °C to 800 °C

The mixtures contain waste tire and *Acrocomia aculeata* endocarp. The filename value `wt` is the waste-tire percentage on a dry-mass basis; the complementary percentage is *A. aculeata* endocarp. Before mixture preparation, the material was oven-dried at 105 °C for at least 48 hours.

No replicate conditions occur in the retained 20-file design. The `2` suffix in three filenames appears to identify replacement or repeated runs, but the corresponding first runs are not present and the reason for replacement is not documented.

### 2.3 Instrument metadata common to the exports

| Field | Confirmed value |
|---|---|
| Export format | NETZSCH5 ANSI, semicolon-separated, decimal point |
| Measurement type | DSC with TG mass channel |
| Instrument | NETZSCH STA 449F3 |
| Project | `copyrolysis` |
| Laboratory | `UNA - FCQ - AI` |
| Operator | `AD` |
| Crucible | DSC/TG pan, Al2O3 |
| Atmosphere | Nitrogen in both purge and protective channels; supplied by Linde Paraguay |
| Purge-channel setting | 80 ml/min in the data columns |
| Protective-channel setting | 20 ml/min in the data columns |
| Temperature calibration | heating-rate-specific 2021 calibration file |
| Correction file | heating-rate-specific correction file |

The exact nitrogen product grade/purity is not recoverable from the exports. The author recalls a grade designation resembling `406` or `408`; this must be reported as unknown unless a cylinder certificate or laboratory record is found.

### 2.4 Raw dynamic columns

| Position | Instrument heading | Meaning | Raw unit | Julia boundary name | Status |
|---:|---|---|---|---|---|
| 1 | `Temp./°C` | measured sample/program temperature | °C | `temperature_C` | Confirmed |
| 2 | `Time/min` | elapsed instrument-program time | min | `time_min` | Confirmed |
| 3 | `DSC/(mW/mg)` | mass-normalized differential scanning calorimetry signal | mW/mg | `dsc_mW_per_mg` | Confirmed |
| 4 | `Mass/%` | instrument-reported remaining mass | % | `mass_percent` | Confirmed |
| 5 | `Gas Flow(purge2)/(ml/min)` | purge-channel volumetric flow setting | ml/min | `purge_flow_mL_per_min` | Confirmed |
| 6 | `Gas Flow(protective)/(ml/min)` | protective-channel volumetric flow setting | ml/min | `protective_flow_mL_per_min` | Confirmed |
| 7 | `Sensit./(uV/mW)` | DSC sensitivity/calibration channel | µV/mW | `sensitivity_uV_per_mW` | Confirmed |
| 8 | `Segment` | instrument program segment identifier | dimensionless integer | `segment` | Confirmed where exported |

The 0 wt% runs at 5, 10, and 15 K/min export seven columns and omit `Segment`; their headers indicate that only segment 2 of 2 was exported. The other 17 runs include segment identifiers 1 and 2.

### 2.5 Run and alias map

The measured rate is the ordinary-least-squares slope of temperature versus time over 100–700 °C. It is an audit diagnostic, not yet the rate that will be used in kinetic calculations.

| Canonical MAT variable | Legacy alias | Raw export | Waste tire (wt%) | Nominal rate (K/min) | Measured rate (K/min) | Sample mass (mg) | Acquisition date |
|---|---|---|---:|---:|---:|---:|---|
| `ExpDat0wt05hr` | `St01` | `ExpDat_0wt05hr.txt` | 0 | 5 | 4.9939 | 10.4455 | 2021-03-18 09:11 |
| `ExpDat25wt05hr2` | `St02` | `ExpDat_25wt05hr2.txt` | 25 | 5 | 4.9938 | 10.4802 | 2021-03-20 03:50 |
| `ExpDat50wt05hr` | `St03` | `ExpDat_50wt05hr.txt` | 50 | 5 | 4.9968 | 10.2436 | 2021-03-17 08:56 |
| `ExpDat75wt05hr` | `St04` | `ExpDat_75wt05hr.txt` | 75 | 5 | 4.9946 | 10.9100 | 2021-03-17 21:05 |
| `ExpDat100wt05hr2` | `St05` | `ExpDat_100wt05hr2.txt` | 100 | 5 | 4.9944 | 10.0911 | 2021-03-19 15:48 |
| `ExpDat0wt10hr` | `St06` | `ExpDat_0wt10hr.txt` | 0 | 10 | 10.0179 | 10.4280 | 2021-03-18 13:55 |
| `ExpDat25wt10hr` | `St07` | `ExpDat_25wt10hr.txt` | 25 | 10 | 10.0238 | 9.9959 | 2021-03-17 00:51 |
| `ExpDat50wt10hr` | `St08` | `ExpDat_50wt10hr.txt` | 50 | 10 | 10.0256 | 10.4044 | 2021-03-17 13:16 |
| `ExpDat75wt10hr` | `St09` | `ExpDat_75wt10hr.txt` | 75 | 10 | 10.0200 | 10.7570 | 2021-03-18 01:20 |
| `ExpDat100wt10hr2` | `St10` | `ExpDat_100wt10hr2.txt` | 100 | 10 | 10.0183 | 10.0453 | 2021-03-19 20:03 |
| `ExpDat0wt15hr` | `St11` | `ExpDat_0wt15hr.txt` | 0 | 15 | 15.1457 | 10.5555 | 2021-03-18 16:53 |
| `ExpDat25wt15hr` | `St12` | `ExpDat_25wt15hr.txt` | 25 | 15 | 15.1715 | 10.1916 | 2021-03-17 03:54 |
| `ExpDat50wt15hr` | `St13` | `ExpDat_50wt15hr.txt` | 50 | 15 | 15.1609 | 10.3115 | 2021-03-17 16:14 |
| `ExpDat75wt15hr` | `St14` | `ExpDat_75wt15hr.txt` | 75 | 15 | 15.1495 | 10.1680 | 2021-03-18 04:18 |
| `ExpDat100wt15hr` | `St15` | `ExpDat_100wt15hr.txt` | 100 | 15 | 15.1463 | 10.0064 | 2021-03-19 23:00 |
| `ExpDat0wt20hr` | `St16` | `ExpDat_0wt20hr.txt` | 0 | 20 | 20.4070 | 10.8195 | 2021-03-18 19:25 |
| `ExpDat25wt20hr` | `St17` | `ExpDat_25wt20hr.txt` | 25 | 20 | 20.4663 | 10.6245 | 2021-03-17 06:31 |
| `ExpDat50wt20hr` | `St18` | `ExpDat_50wt20hr.txt` | 50 | 20 | 20.4301 | 10.3153 | 2021-03-17 18:46 |
| `ExpDat75wt20hr` | `St19` | `ExpDat_75wt20hr.txt` | 75 | 20 | 20.4069 | 10.3628 | 2021-03-18 06:51 |
| `ExpDat100wt20hr` | `St20` | `ExpDat_100wt20hr.txt` | 100 | 20 | 20.3984 | 10.3815 | 2021-03-20 01:31 |

Every mapping above was verified by exact element-by-element comparison, including `NaN` positions, between `code/2021corridas/2021corridas.mat` and `code/2021mar.mat`.

### 2.6 Canonical Julia boundary representation

Import must preserve raw columns and metadata. Numerical routines will receive explicitly converted fields:

```text
temperature_K = temperature_C + 273.15
time_min      = raw time in minutes
mass_fraction = mass_percent / 100
waste_tire_mass_fraction = declared_wt_percent / 100
acrocomia_endocarp_mass_fraction = 1 - waste_tire_mass_fraction
```

Whether `time_min` should be re-zeroed at the ramp start is a preprocessing decision. It is not part of raw import.

## 3. Legacy derived dynamic datasets

### `code/2021marcorrected.mat`

Contains cropped four-column arrays. Values indicate `[temperature_C, time_min, DSC_or_rate, mass_percent]`, but the transformations are not fully documented. Treat as a comparison artifact.

### `code/corridasjulio.mat` and `code/julio2.mat`

Contain five-column `St01`–`St20` arrays used by much of the final analysis code. Their apparent layout is `[temperature_C, time_min, mass_percent, sensitivity_or_placeholder, segment]`. These arrays are normalized or otherwise transformed relative to the raw exports. In the 20 K/min block, column 4 is a constant placeholder. They are not suitable as primary raw data.

### `code/a*neum.mat` and `code/*_nuevo.mat`

Saved MATLAB workspaces containing raw-like arrays mixed with fitted parameters, processed curves, function handles, residuals, and plot inputs. They are result snapshots only.

## 4. Ramp-and-hold validation dataset

### 4.1 Selected validation dataset

`code/isotherm.mat` contains 15 arrays named by composition and nominal hold temperature:

```text
b{0,25,50,75,100}t{450,500,550}
```

This is a complete 5-composition × 3-hold-temperature design. These are not strictly isothermal experiments: each program contains conditioning, a dynamic ramp, and a final temperature hold. Their intended purpose was to test whether kinetic parameters estimated from the dynamic 5/10/15/20 K/min runs can predict a different temperature program. They are therefore reserved for validation and must not be used to fit those parameters.

The 450/500/550 labels are nominal hold temperatures in °C. Numerical inspection of all 15 runs shows:

- segment 3 ramps from approximately 100 °C toward the hold temperature at 20.6877–20.9008 K/min
- segment 4 holds the target for 59.90–59.95 minutes
- the final 200 samples average within 0.023–0.044 °C above the nominal target

`simuisotherm.m`, `movietga.m`, and `somethingsomething.m` use these arrays.

The eight-column layout is strongly consistent with:

| Position | Provisional meaning | Provisional unit | Confidence |
|---:|---|---|---|
| 1 | temperature | °C | High |
| 2 | time | min | High |
| 3 | DSC signal | probably mW/mg | Medium |
| 4 | remaining mass | % | High; corroborated by `cleandatami.m` |
| 5 | purge-channel flow | ml/min | Medium |
| 6 | protective-channel flow | ml/min | Medium |
| 7 | sensitivity | probably µV/mW | Medium |
| 8 | segment identifier, values 1–4 | dimensionless integer | High |

Original instrument text exports or headers have not been located, so the channel names in positions 3 and 5–7 remain provisional. Temperature, time, mass, segment identity, target temperatures, and the validation role are sufficiently established for planning.

### 4.2 Alternative isothermal arrays

`code/2021apr.mat` and `code/2021apr2.mat` contain nine arrays named `m{25,50,75}isot{450,500,550}`. They use a different nine-column corrected-data layout and only cover the three mixtures. These are not selected for predictive validation. Their relationship to `isotherm.mat` may be investigated later as legacy provenance work.

## 5. Conversion and rate quantities

### 5.1 Raw mass

The instrument reports remaining mass as a percentage. Mixture preparation used dry-mass fractions after oven conditioning at 105 °C for at least 48 hours. The instrument's internal normalization still requires care: the first recorded point is usually near 100%, but the 0 wt% exports that begin at segment 2 start around 94%.

### 5.2 Legacy dynamic conversion

The most-used preprocessing function, `cleandata.m`, retained 150–700 °C and defined

```text
alpha = (m_initial - m) / (m_initial - m_final)
```

where `m_initial` and `m_final` were the first and last retained samples.

The accepted primary convention for the rewrite is a thesis-compatible 150–700 °C conversion interval:

```text
alpha_700(T) = [m(150 °C) - m(T)] / [m(150 °C) - m(700 °C)]
```

The reference masses will be evaluated at the exact temperatures using a tested local robust fit or shape-preserving interpolation, rather than whichever rows happen to be first and last. The preprocessing report must also repeat the analysis with alternative final references at 650 °C, 750 °C, and the highest common observed temperature without extrapolation. A 120 °C versus 150 °C initial-reference sensitivity will quantify any residual drying-region effect.

`alpha_700` is a declared analysis coordinate, not a claim that all chemical reactions cease at 700 °C. Raw mass and residual mass above 700 °C remain available for reporting.

### 5.3 Rate definitions

The rewrite will distinguish:

- `dalpha_dt`, with unit min⁻¹
- `dalpha_dT`, with unit K⁻¹
- `dmass_dt`, with mass-basis unit per minute
- `dmass_dT`, with mass-basis unit per kelvin

For an ideal linear ramp with heating rate `beta = dT/dt`:

```text
dalpha_dt = beta * dalpha_dT
```

Legacy four-column processed arrays do not use column 3 consistently; its meaning must never be inferred from matrix width alone.

## 6. Resolved decisions and remaining provenance questions

| ID | Resolved decision | Evidence |
|---|---|---|
| Q001 | The non-tire component is *Acrocomia aculeata* endocarp. | Author confirmation, 2026-08-17 |
| Q002 | Mixtures are expressed as dry-mass percentages after oven drying at 105 °C for at least 48 hours. | Author confirmation, 2026-08-17 |
| Q003 | Purge and protective atmospheres used nitrogen supplied by Linde Paraguay. | Author confirmation, 2026-08-17 |
| Q004 | Primary conversion uses the 150–700 °C interval, with endpoint-sensitivity analysis. | Original intent plus accepted rewrite decision, 2026-08-17 |
| Q008 | `code/isotherm.mat` is the ramp-and-hold predictive-validation set; its labels are hold temperatures in °C. | Author confirmation and numerical program audit, 2026-08-17 |
| Q010 | Scientific acceptance requires out-of-sample prediction overlays and residuals; no legacy figure is a pixel-matching target. | Rewrite validation objective, 2026-08-17 |

Non-blocking provenance questions retained for later investigation:

| ID | Question | Treatment until resolved |
|---|---|---|
| Q003a | What exact Linde nitrogen grade or purity was used? | Record nitrogen and supplier; mark exact grade unknown |
| Q005 | Which buoyancy/baseline corrections were applied inside the NETZSCH export? | Preserve raw channels; compare correction strategies during M3 |
| Q006 | Why were three dynamic filenames suffixed `2`? | Treat retained files as selected runs, not replicate pairs |
| Q007 | Are unretained first runs or additional replicates available? | Do not claim experimental replicate uncertainty |
| Q009 | Why do ramp-and-hold flow columns contain 20/20 rather than dynamic-run 80/20? | Report recorded values and test cross-program sensitivity if needed |

## 7. M0 data gate

The M0 data gate is satisfied for beginning the Julia foundation and read-only importer. Dynamic source columns, units, sample identities, composition basis, run aliases, heating rates, conversion convention, and validation-program role are documented. Remaining questions concern historical provenance and must stay visible, but they do not block M1 or M2.

## 8. M2 implementation status

`config/datasets.toml` now maps all 20 dynamic and 15 ramp-and-hold experiments to explicit
MAT variables and fingerprinted source files. The Julia importer retains all raw matrix rows,
converts the temperature boundary to kelvin, validates structure and metadata, and records
the known trailing invalid rows without removing them. Dynamic entries also retain their
authoritative NETZSCH text-export paths.

The reproducible audit is available in `docs/data_audit.md`; its machine-readable counterpart
is `docs/audits/m2_data_inventory.toml`. Row filtering, time re-zeroing, segment selection,
mass-reference estimation, conversion, smoothing, and derivatives remain M3 responsibilities.
