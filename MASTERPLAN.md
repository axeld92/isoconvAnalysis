# Isoconversional Analysis Julia Rewrite — Masterplan

Last updated: 2026-08-17  
Project status: M6 complete — M7 ready to start  
Target language: Julia 1.11  
Legacy implementation: `code/` (read-only reference)

## 1. Purpose

Build a trustworthy, reproducible Julia implementation of the thesis workflow while preserving its useful scientific ideas:

1. Import and characterize thermogravimetric analysis (TGA) experiments.
2. Convert mass-loss measurements into conversion and conversion-rate curves.
3. Estimate conversion-dependent activation energy with multiple isoconversional methods.
4. Deconvolve overlapping DTG peaks with Fraser–Suzuki functions.
5. Estimate kinetic parameters and candidate reaction models.
6. Simulate non-isothermal and isothermal conversion histories.
7. Quantify mixture interaction or synergy.
8. Generate publication-quality diagnostics and reproducible reports.

The rewrite will preserve the research questions, not undocumented implementation behavior. Legacy numerical results are comparison artifacts, not ground truth.

## 2. Success criteria

The rewrite is successful when:

- [ ] A new user can instantiate the Julia environment and reproduce an analysis from one documented command.
- [ ] Every input column, internal variable, and reported quantity has an explicit meaning and unit.
- [ ] Raw measurements remain immutable and can always be traced to their source file and variable.
- [ ] Preprocessing decisions are configurable, visible in QC plots, and covered by sensitivity analysis.
- [x] Each isoconversional method passes synthetic tests with known activation energy.
- [x] Advanced Vyazovkin calculations work for any supported number of experiments without hard-coded objective constants.
- [x] Confidence intervals have a documented statistical meaning and are never replaced by optimizer residuals.
- [x] Deconvolution results are invariant to experiment input order within numerical tolerance.
- [x] Peak count is supported by information criteria and structural diagnostics rather than manual relabeling.
- [ ] Cross-rate component identity is validated well enough for kinetic propagation.
- [ ] Simulations enforce physically valid conversion states instead of discarding imaginary parts.
- [ ] At least one heating program is withheld from fitting and used for predictive validation.
- [ ] Tests, documentation, examples, and a locked Julia environment are included.

## 3. Guiding principles

1. **Units are part of the API.** Temperature, time, mass, heating rate, activation energy, and rate variables must be unambiguous.
2. **Raw data are immutable.** Corrections and transformations create new data products.
3. **Scientific functions are pure.** Computation does not implicitly load workspace variables, create figures, or mutate inputs.
4. **Plotting is downstream.** Numerical routines return structured results and diagnostics; plotting functions consume them.
5. **Validation precedes feature breadth.** A smaller validated implementation is preferable to a complete unverified port.
6. **Failures remain visible.** Invalid domains, non-monotonic conversion, non-identifiable fits, and optimizer failures produce explicit diagnostics.
7. **Deconvolution is descriptive until validated.** A fitted peak is not automatically a physical reaction or chemical component.
8. **Legacy behavior is opt-in.** Any compatibility mode must be named, documented, and tested separately.

## 4. Scope

### In scope for version 1.0

- Dynamic TGA runs at multiple heating rates.
- Existing isothermal datasets, after their schema is confirmed.
- Friedman, KAS, FWO, Starink, and advanced Vyazovkin analyses.
- Fraser–Suzuki peak models with model selection and uncertainty diagnostics.
- Candidate Sesták–Berggren reaction models.
- Non-isothermal and ramp-plus-hold simulations.
- Mixture-additivity and synergy analysis.
- Static reports and publication-quality figures.
- Julia package, command-line scripts, tests, and documentation.

### Initially out of scope

- A graphical user interface.
- Automatic interpretation of unnamed MAT-file columns.
- Exact replication of every legacy figure or saved workspace value.
- Porting unrelated experiments, tutorial scripts, or the vendored TGAnalysis project.
- Treating peak labels such as “cellulose” or “lignin” as established without independent evidence.
- Reactor-scale heat and mass transfer modeling.

## 5. Repository strategy

The existing `code/` directory remains untouched as a legacy archive. The Julia implementation will live at the repository root.

```text
ISOCONVERSIONAL_ANALYSIS/
├── Project.toml
├── Manifest.toml
├── README.md
├── MASTERPLAN.md
├── src/
│   ├── IsoconversionalAnalysis.jl
│   ├── Types.jl
│   ├── DataIO.jl
│   ├── Preprocessing.jl
│   ├── Isoconversional.jl
│   ├── PeakModels.jl
│   ├── Deconvolution.jl
│   ├── ReactionModels.jl
│   ├── Simulation.jl
│   ├── Synergy.jl
│   ├── Validation.jl
│   ├── Plotting.jl
│   └── Reporting.jl
├── test/
│   ├── runtests.jl
│   ├── fixtures/
│   ├── test_data_io.jl
│   ├── test_preprocessing.jl
│   ├── test_isoconversional.jl
│   ├── test_deconvolution.jl
│   ├── test_reaction_models.jl
│   └── test_simulation.jl
├── scripts/
│   ├── audit_legacy_data.jl
│   ├── run_analysis.jl
│   └── reproduce_reference_report.jl
├── config/
│   ├── datasets.toml
│   └── analysis_defaults.toml
├── data/
│   ├── README.md
│   ├── interim/
│   └── processed/
├── results/                 # generated and ignored by Git
├── docs/
│   ├── data_dictionary.md
│   ├── methods.md
│   ├── validation.md
│   └── legacy_mapping.md
└── code/                    # existing MATLAB archive; never imported implicitly
```

Large raw MAT files may remain in `code/` initially. `config/datasets.toml` will map stable experiment identifiers to source files and MAT variables without copying or modifying the originals.

## 6. Core data contracts

### 6.1 Raw experiment

Define an `Experiment` type containing, at minimum:

- `id::String`
- `source_file::String`
- `source_variable::String`
- `sample_id::String`
- `composition` with named component fractions
- `temperature_K::Vector{Float64}`
- `time_min::Vector{Float64}`
- `mass` with a documented raw unit
- `nominal_heating_rate_K_per_min::Float64`
- experimental metadata: atmosphere, flow rate, initial mass, date, replicate, notes
- import warnings and provenance metadata

Arrays used in numerical solvers will use canonical units named in their fields. Import and reporting boundaries must perform explicit conversions. No function may infer Celsius versus Kelvin from numeric magnitude.

### 6.2 Processed experiment

Define a `ProcessedExperiment` containing:

- a reference to the raw experiment identity
- retained temperature and time vectors
- normalized mass fraction
- conversion `alpha`
- `dalpha_dT` and, when needed, `dalpha_dt`
- measured or fitted heating-rate curve
- preprocessing configuration
- QC metrics and warnings
- hashes or identifiers connecting the result to source data and configuration

### 6.3 Analysis results

Use dedicated result types rather than anonymous matrices:

- `IsoconversionalResult`
- `DeconvolutionResult`
- `ReactionModelResult`
- `SimulationResult`
- `SynergyResult`

Each result must include estimates, uncertainty representation, diagnostics, method name, configuration, convergence status, and provenance.

## 7. Scientific specifications

### 7.1 Data ingestion and audit

- Read MATLAB v5 files with explicit variable selection.
- Build a data dictionary before assigning semantic meaning to each column.
- Verify finite values, vector lengths, time ordering, temperature ordering, sampling intervals, and segment markers.
- Detect trailing NaN rows without assuming there is only one.
- Estimate the actual heating-rate curve from temperature versus time and compare it with the nominal rate.
- Record rather than silently remove duplicated, reversed, or non-monotonic points.
- Export a machine-readable audit summary and a human-readable data report.

Required gate: the meanings of the columns in the selected `St01`–`St20` source dataset and the mapping from run number to composition/heating rate must be confirmed before kinetic calculations begin.

### 7.2 Preprocessing

The preprocessing configuration must explicitly define:

- temperature or time interval selection
- baseline and buoyancy correction, if applicable
- how initial and final masses are estimated
- smoothing method and its window in physical or clearly documented sampling units
- derivative method
- monotonic-conversion policy
- valid conversion range for later analysis

Default workflow to evaluate:

1. Remove invalid acquisition rows.
2. Select the analysis segment using recorded segment metadata where available.
3. Estimate initial and final mass from stable windows rather than single boundary points.
4. Calculate conversion as
   `alpha = (m_initial - m) / (m_initial - m_final)`.
5. Diagnose values outside `[0, 1]`; do not silently clip them.
6. Smooth mass or conversion only when justified by a sensitivity study.
7. Differentiate with a local polynomial/Savitzky–Golay or other validated method.
8. Verify that integrating `dalpha_dT` reconstructs conversion within tolerance.
9. Generate a QC figure showing raw mass, corrected mass, conversion, derivative, retained interval, and residual reconstruction error.

### 7.3 Isoconversional analysis

Implement all methods through one consistent interface:

```julia
analyze(method, experiments, alpha_grid; options...)
```

Every method must:

- use absolute temperature internally
- interpolate each experiment at the same conversion values
- reject or flag extrapolated conversion values
- support an arbitrary valid number of heating programs where the method permits it
- return per-conversion diagnostics and uncertainty metadata
- avoid creating figures or relying on global variables

Methods:

- **Friedman:** differential method using `log(dalpha/dt)` or the equivalent consistently derived rate.
- **KAS:** linear integral approximation.
- **FWO:** linear integral approximation with its documented coefficient and logarithm convention.
- **Starink:** linear integral approximation with documented exponent and coefficient.
- **Advanced Vyazovkin:** incremental numerical integration over a configurable conversion interval.

Advanced Vyazovkin requirements:

- Compute the pairwise objective from the actual number of experiments; never hard-code `12`.
- Minimize the objective itself rather than distance from an assumed theoretical value.
- Use a bounded one-dimensional optimizer over a configurable physical energy interval.
- Use stable numerical integration and guard against underflow/overflow.
- Make `delta_alpha` and integration resolution explicit.
- Implement a literature-backed confidence-interval procedure separately from optimization quality metrics.
- Report boundary solutions and flat or multimodal objectives.

### 7.4 Uncertainty

Uncertainty outputs must state exactly what they represent. Planned layers:

1. Regression coefficient intervals for linear isoconversional methods.
2. Method-specific intervals for the advanced Vyazovkin method.
3. Sensitivity to preprocessing configuration and conversion interval.
4. Parametric or residual bootstrap where statistically defensible.
5. Across-replicate variation if replicate experiments are available.

Optimizer residuals, objective values, and fit MSE are diagnostics, not confidence intervals.

### 7.5 Fraser–Suzuki deconvolution

Implement one numerically stable Fraser–Suzuki kernel with a Gaussian limit as skew approaches zero.

Requirements:

- positive peak heights and widths
- ordered peak centers to prevent label switching
- explicit domain handling for the logarithm
- configurable two-, three-, and four-peak models
- deterministic fitting for a fixed seed and configuration
- invariance to the order in which heating-rate experiments are supplied
- residual plots and integrated-area reconstruction checks
- parameter covariance or bootstrap intervals where identifiable
- AICc/BIC and residual diagnostics for peak-count comparison

Evaluate two strategies:

1. Independent fits at each heating rate followed by component matching.
2. Joint fitting across heating rates with shared component identity and selected shared or rate-dependent parameters.

The joint approach is preferred if synthetic and residual tests show that it is identifiable.

### 7.6 Kinetic triplet and compensation analysis

Preserve the original compensation-effect and Sesták–Berggren ideas, but treat them as a later, separately validated layer.

- Estimate pre-exponential factors without duplicating non-identifiable intercept terms.
- Quantify the correlation between `E` and `log(A)`.
- Fit candidate reaction functions only on a documented conversion interval.
- Compare candidate models using prediction error and information criteria.
- Keep empirical model identification distinct from mechanistic claims.
- Propagate uncertainty from activation energy and deconvolution where feasible.

### 7.7 Simulation and prediction

- Implement dynamic temperature programs as explicit functions `T(t)`.
- Support linear ramps, isothermal holds, and ramp-plus-hold programs.
- Integrate with DifferentialEquations.jl or an equivalent validated solver.
- Enforce physical state bounds through formulation, callbacks, or termination conditions.
- Never use `real(...)` to conceal invalid model states.
- Compare predictions with held-out dynamic and isothermal experiments.
- Report residuals in conversion and rate, not only a single aggregate MSE.

### 7.8 Synergy analysis

- Define the additive reference using explicit mass or composition fractions.
- Interpolate component experiments onto a common physical domain without uncontrolled extrapolation.
- Report absolute and relative deviations in mass, conversion, and conversion rate as distinct quantities.
- Include integrated deviation, peak deviation, and uncertainty bands.
- Document whether a positive value means acceleration, additional conversion, or another effect.

## 8. Validation strategy

### Level 1 — Unit tests

- Known values for unit conversion and conversion normalization.
- Analytical derivatives on simple functions.
- Fraser–Suzuki Gaussian-limit and domain tests.
- Objective-function permutation invariance.
- Serialization and configuration round trips.

### Level 2 — Synthetic kinetic tests

- Single-step first-order process with constant known activation energy.
- Sesták–Berggren cases with known parameters.
- Multiple heating rates with exact and noisy data.
- Conversion-dependent activation-energy cases.
- Two- and three-peak synthetic mixtures with known areas and centers.
- Irregular sampling, missing points, and mild temperature-program deviations.

Acceptance targets will be stated before generating each fixture so they cannot be tuned after seeing results.

### Level 3 — Legacy-data regression tests

- Select a small, versioned subset of original runs as fixtures.
- Verify import, preprocessing invariants, and deterministic output.
- Compare corrected Julia outputs with legacy outputs while explaining expected differences.
- Never change corrected code merely to reproduce a known legacy defect.

### Level 4 — Predictive validation

- Fit with a subset of heating rates and predict a held-out rate.
- Fit dynamic data and predict the isothermal experiments.
- Compare conversion curves, rate curves, peak locations, residual structure, and uncertainty coverage.

## 9. Milestones and progress tracker

Status legend: `Not started` · `In progress` · `Blocked` · `Complete`

| ID | Milestone | Status | Exit evidence |
|---|---|---|---|
| M0 | Scientific requirements and data dictionary | Complete | Confirmed dataset map and documented equations |
| M1 | Julia package scaffold and reproducible environment | Complete | Package loads; 31 tests and local CI-equivalent gates pass |
| M2 | MAT ingestion and legacy data audit | Complete | All 35 runs load; deterministic Markdown/TOML audit and 89 tests pass |
| M3 | Validated preprocessing pipeline | Complete | QC report, sensitivity artifacts, and 158 passing tests |
| M4 | Isoconversional methods | Complete | Synthetic gates, real-data audit, figures, and 251 tests pass |
| M5 | Deconvolution framework | Complete | Synthetic recovery, order invariance, real-data audit, and 310 tests pass |
| M6 | Kinetic triplet estimation | Complete | Synthetic recovery, model comparison, real-data predictive audit, and 369 tests pass |
| M7 | Dynamic/ramp-and-hold simulation | Not started | Held-out prediction report |
| M8 | Synergy analysis and publication reporting | Not started | Reproducible reference report |
| M9 | Legacy comparison, documentation, and v1.0 release | Not started | Tagged release and archived validation bundle |

### M0 — Scientific requirements and data dictionary

- [x] Select the canonical raw MAT dataset. See [`docs/data_dictionary.md`](docs/data_dictionary.md).
- [x] Confirm the meaning and unit of every source column in the canonical dynamic dataset.
- [x] Confirm run-to-composition and run-to-heating-rate mapping by exact array comparison.
- [x] Identify retained dynamic runs as distinct conditions; no replicate pair is present in the canonical set.
- [x] Confirm dry-mass mixture basis and adopt a primary 150–700 °C conversion coordinate with sensitivity checks.
- [x] Transcribe the recovered equations with primary references. See [`docs/methods.md`](docs/methods.md).
- [x] Record expected scientific outputs; ramp-and-hold prediction overlays and residuals are mandatory validation evidence.
- [x] Classify legacy scripts as core, reference, exploratory, or unrelated. See [`docs/legacy_mapping.md`](docs/legacy_mapping.md).

M0 conclusion: the feedstocks, dry-mass basis, nitrogen atmosphere, conversion convention, dynamic experiment map, and ramp-and-hold validation role are documented. Exact gas grade and historical correction provenance remain visible as non-blocking metadata gaps.

Exit gate: no numerical implementation starts with an unresolved source-data column or unit.

### M1 — Julia foundation

- [x] Initialize `Project.toml` for Julia 1.11 with locked runtime, test, and tooling environments.
- [x] Create the package/module structure.
- [x] Add formatting, linting, documentation, and test commands.
- [x] Add `.gitignore` rules for generated results and large local artifacts.
- [x] Establish configuration parsing and structured logging.
- [x] Add continuous integration for supported Julia versions.
- [x] Record dependency choices and licenses in [`docs/dependencies.md`](docs/dependencies.md).

M1 conclusion: the package loads on Julia 1.11.6, the default M0 configuration validates,
31 tests including Aqua pass, JuliaFormatter reports a clean tree, and Documenter builds the
manual without errors. Runtime and development dependencies are isolated and pinned. See
[`docs/m1_foundation.md`](docs/m1_foundation.md) for commands and evidence. Remote CI execution
will become observable once the repository has a configured remote.

Exit gate: a fresh clone can instantiate, load the package, and run the baseline test suite.

### M2 — Data ingestion

- [x] Implement explicit MAT-file and variable loading with source checksum verification.
- [x] Implement `Experiment` validation while preserving and reporting invalid raw rows.
- [x] Create [`config/datasets.toml`](config/datasets.toml) for 20 calibration and 15 validation runs.
- [x] Generate a complete machine-readable inventory in [`docs/audits/m2_data_inventory.toml`](docs/audits/m2_data_inventory.toml).
- [x] Add small representative MAT fixtures and real-source integration tests.
- [x] Produce the first human-readable [`docs/data_audit.md`](docs/data_audit.md) report.

M2 conclusion: all 35 selected MAT variables load without workspace preparation, both source
hashes match, no configured variable is missing, and all experiments pass structural
validation. The importer retains the 35 known trailing invalid rows and exposes them through
warnings rather than silently preprocessing them. The audit command produces byte-identical
Markdown and TOML outputs across reruns. The complete suite passes 89 tests, including Aqua.

Exit gate: all selected raw experiments load without manual workspace preparation.

### M3 — Preprocessing

- [x] Implement interval selection and invalid-row handling.
- [x] Implement configurable mass-baseline estimation.
- [x] Implement conversion calculation.
- [x] Compare candidate smoothing and derivative methods.
- [x] Implement integration-based derivative consistency checks.
- [x] Add preprocessing QC plots and warnings.
- [x] Freeze a default configuration only after sensitivity review.

M3 conclusion: raw experiments now transform into source-row-traceable `ProcessedExperiment`
records without mutation, clipping, or silent monotonic repair. Synthetic polynomial and
end-to-end tests validate differentiation and integration. The complete real-data
[`preprocessing audit`](docs/preprocessing_audit.md) compares derivative, smoothing,
reference-mass, and endpoint choices across the calibration design; the accepted ±5 K cubic
local-polynomial profile keeps the maximum reconstruction RMSE below `1.22e-4`. CairoMakie QC
figures remain downstream in the documentation environment, and all ramp-and-hold programs
remain held out.

Exit gate: preprocessing passes synthetic tests and produces an approved QC report on real data.

### M4 — Isoconversional methods

- [x] Implement shared interpolation and regression utilities.
- [x] Implement Friedman.
- [x] Implement KAS.
- [x] Implement FWO.
- [x] Implement Starink.
- [x] Implement advanced Vyazovkin.
- [x] Implement uncertainty and diagnostics.
- [x] Cross-check methods on synthetic and real data.
- [x] Document formula conventions and references.

M4 conclusion: one typed, order-invariant API now implements all five methods on a common
first-upward-crossing conversion grid. Known-energy and noisy first-order fixtures satisfy
predeclared method-specific tolerances. Advanced Vyazovkin minimizes the true pairwise
objective for arbitrary experiment counts with log-scaled quadrature, bounded optimization,
and explicit boundary/flat/multimodal diagnostics. Linear methods return transformed
Student-t slope intervals; the advanced method returns Vyazovkin–Wight pairwise-variance
Fisher intervals. The complete real-data [`isoconversional audit`](docs/isoconversional_audit.md)
reports all 25 method/composition profiles and 100 endpoint-sensitivity comparisons without
using the held-out ramp-and-hold experiments. The full package suite passes 251 tests,
including Aqua.

Exit gate: known-energy synthetic cases meet predeclared accuracy tolerances.

### M5 — Deconvolution

- [x] Implement and test the Fraser–Suzuki kernel.
- [x] Implement constrained single-experiment fits.
- [x] Implement experiment-order invariance tests.
- [x] Implement peak-count comparison.
- [x] Evaluate joint multi-rate fitting.
- [x] Implement component-area and reconstruction checks.
- [x] Add parameter uncertainty and identifiability diagnostics.

M5 conclusion: a stable Fraser–Suzuki kernel and typed independent/joint result contracts now
fit measured `dalpha/dT` without clipping. Smooth transformations enforce nonnegative heights,
bounded shape and width, and ordered separated centers. Known three-peak fixtures meet the
predeclared center, area-fraction, and RMSE tolerances; noisy fixtures select the correct peak
count; and joint fits are invariant to caller order. The complete real-data
[`deconvolution audit`](docs/deconvolution_audit.md) retains all 60 two-/three-/four-peak
candidates for the 20 calibration runs. Raw BIC minima are three peaks in 8 runs and four in
12; final selections after structural filtering are two peaks in 4 runs, three in 10, and four
in 6. One run has no fully eligible candidate and retains its raw minimum with an explicit
warning. The evaluated
shared-skew joint model fails the convergence or 5% objective-loss gate for some compositions,
so independent fits remain the M5 default and cross-rate component identity is not yet assumed.
Audit Markdown, TOML, and SVG artifacts are byte-identical across reruns. The full package
suite passes 310 tests, including Aqua; JuliaFormatter and Documenter gates pass.

Exit gate: synthetic peaks and areas are recovered within predeclared tolerances, with failures reported clearly.

### M6 — Kinetic triplet

- [x] Implement 17 fixed and two empirical candidate reaction functions in one named registry.
- [x] Estimate the pre-exponential factor in `min^-1` from one identifiable log-rate intercept.
- [x] Implement two- and three-exponent Šesták–Berggren fits without a rank-deficient
  parameterization.
- [x] Evaluate the compensation-effect assumptions without substituting same-data correlation
  for predictive evidence.
- [x] Add BIC/AICc comparison, structural filtering, conditional Student-t uncertainty, and
  leave-one-heating-rate-out validation.

M6 conclusion: the overall-conversion API recovers known `(E, ln A, m, n)` synthetic
parameters within `1e-5`, selects the expected noisy empirical model, exposes nonpositive
rates, and is invariant to experiment order. The complete real-data
[`reaction-model audit`](docs/reaction_model_audit.md) retains 95 candidate fits, 380
heating-rate validation rows, and 20 single-rate compensation diagnostics. Raw criterion
minima, structural screening, covariance/correlation, residual autocorrelation, and
compensation extrapolation remain visible. Every selected composition-level constant-triplet
model fails the predeclared held-out-rate log-RMSE gate, so these fits are descriptive M7
baselines rather than validated predictive mechanisms. No M5 peak is assigned a common
reaction identity. Audit Markdown, TOML, and SVG artifacts are byte-identical across reruns.
The full package suite passes 369 tests, including Aqua; JuliaFormatter and Documenter gates
pass.

Exit gate: parameters are recoverable from synthetic cases, and every real-data candidate is
either predictively accepted or explicitly rejected before simulation.

### M7 — Simulation

- [ ] Implement temperature-program types.
- [ ] Implement non-isothermal ODE prediction.
- [ ] Implement isothermal and ramp-plus-hold prediction.
- [ ] Add state-domain enforcement and solver diagnostics.
- [ ] Perform held-out heating-rate validation.
- [ ] Perform isothermal validation.

Exit gate: predictions and residuals are reproducible, finite, and physically valid.

### M8 — Synergy and reporting

- [ ] Implement additive mixture references.
- [ ] Implement synergy metrics and uncertainty bands.
- [ ] Select a publication plotting backend.
- [ ] Produce standard QC, activation-energy, deconvolution, and prediction figures.
- [ ] Generate a reproducible HTML/PDF or Markdown reference report.

Exit gate: one command regenerates the complete reference analysis and its provenance record.

### M9 — Release

- [ ] Complete legacy-versus-Julia comparison.
- [ ] Document intentional deviations from the thesis code.
- [ ] Complete user and developer documentation.
- [ ] Review dependency licenses and data-sharing constraints.
- [ ] Run the full validation matrix from a clean environment.
- [ ] Tag version 1.0 and archive the validation outputs.

## 10. Immediate next sprint — M7

- [ ] Define typed linear-ramp, isothermal-hold, and ramp-plus-hold temperature programs.
- [ ] Implement a physical-domain ODE interface with explicit solver and termination
  diagnostics.
- [ ] Retain M6 constant-triplet fits only as frozen descriptive benchmarks.
- [ ] Construct a conversion-dependent `E(alpha)` and `A(alpha)` path from the validated M4
  results without circular use of the held-out data.
- [ ] Verify synthetic dynamic and ramp-plus-hold trajectories before using real data.
- [ ] Predict all 15 `code/isotherm.mat` programs without refitting and report ramp/hold
  residuals separately.
- [ ] Compare constant-triplet and conversion-dependent predictions using predeclared error
  and physical-validity gates.

M7 starts with a negative but useful M6 result: one constant triplet cannot predict the
dynamic heating rates adequately. The simulation layer must expose this benchmark while
testing whether conversion-dependent isoconversional parameters improve genuinely held-out
ramp-plus-hold predictions.

## 11. Proposed Julia dependencies

Dependencies will be added only when their role is established. Initial candidates:

| Purpose | Candidate packages |
|---|---|
| MAT-file input | MAT.jl |
| Tables and export | DataFrames.jl, CSV.jl |
| Configuration | TOML standard library |
| Interpolation | Interpolations.jl or DataInterpolations.jl |
| Optimization | Optim.jl and/or Optimization.jl |
| Nonlinear least squares | LsqFit.jl or a selected Optimization.jl backend |
| Automatic differentiation | ForwardDiff.jl |
| Differential equations | DifferentialEquations.jl |
| Statistics | StatsBase.jl, Distributions.jl, GLM.jl as needed |
| Figures | CairoMakie.jl |
| Testing | Test standard library, Aqua.jl, JET.jl |
| Documentation | Documenter.jl |

Package selection must consider maintenance status, reproducibility, numerical behavior, and transitive dependency cost.

## 12. Decision log

| ID | Decision | Status | Rationale |
|---|---|---|---|
| D001 | Use Julia 1.11 as the primary runtime | Accepted | Julia 1.11.6 is available in the working environment |
| D002 | Preserve `code/` unchanged as the legacy archive | Accepted | Prevent accidental loss and keep historical traceability |
| D003 | Build a Julia package, not a collection of workspace scripts | Accepted | Enables testing, reuse, documentation, and deterministic execution |
| D004 | Use canonical internal units with unit-bearing field names and boundary validation | Accepted | Prevents Celsius/Kelvin and `dalpha/dt`/`dalpha/dT` drift without complicating every solver |
| D005 | Treat saved MATLAB results as comparison artifacts | Accepted | The audit found calculation and uncertainty defects |
| D006 | Validate on synthetic data before real-data interpretation | Accepted | Provides known truth for numerical methods |
| D007 | Keep deconvolution and mechanistic interpretation separate | Accepted | Peak fits alone do not identify chemical mechanisms |
| D008 | Prefer configuration files over hard-coded thresholds | Accepted | Makes decisions reviewable and reproducible |
| D009 | Use a joint multi-rate deconvolution if identifiability tests support it | Superseded by D043 | M5 evaluated the proposal; the real-data gate does not support joint fitting as the default |
| D010 | Use CairoMakie for publication graphics | Accepted | M3 adopts it in the isolated documentation environment for high-quality static QC output |
| D011 | Treat the 20 NETZSCH text exports as the primary dynamic source and `code/2021corridas/2021corridas.mat` as the canonical MAT snapshot | Accepted | The text files retain instrument metadata; the MAT arrays match them and support regression tests |
| D012 | Treat `code/2021mar.mat` as alias evidence and later saved workspaces as derived comparison artifacts | Accepted | Exact comparison recovers the `St01`–`St20` mapping without confusing processed workspaces with raw data |
| D013 | Identify the feedstocks as waste tire and *Acrocomia aculeata* endocarp, mixed by dry mass after at least 48 h at 105 °C | Accepted | Author-confirmed sample preparation |
| D014 | Record nitrogen in both gas channels and Linde Paraguay as supplier; leave exact grade unknown | Accepted | Gas identity is confirmed but the historical grade designation is uncertain |
| D015 | Use 150–700 °C as the primary conversion coordinate and require reference-temperature sensitivity analysis | Accepted | Preserves the thesis intent while exposing dependence on normalization endpoints |
| D016 | Reserve all 15 `code/isotherm.mat` ramp-and-hold runs for out-of-sample validation | Accepted | They were designed to test prediction under temperature programs not used for parameter estimation |
| D017 | Reproduce scientific outputs and diagnostics, not legacy figures pixel-for-pixel | Accepted | Validation quality is more meaningful than matching historical plot styling |
| D018 | Keep the M1 runtime standard-library-only and isolate test and documentation tools in separately locked environments | Accepted | Minimizes the scientific runtime while preserving exact development-tool versions |
| D019 | Resolve relative output paths from the configuration file and reject noncanonical boundary-unit labels | Accepted | Makes runs independent of the shell working directory and prevents silent unit drift |
| D020 | Use JuliaFormatter, Aqua, Documenter, and Julia 1.11 GitHub Actions as mandatory foundation gates | Accepted | Provides deterministic style, package hygiene, checked API documentation, and repeatable automation |
| D021 | Use MAT.jl for checksum-verified, explicitly named variable access to the canonical MATLAB snapshots | Accepted | Satisfies the legacy v5 import requirement without loading arbitrary workspace contents |
| D022 | Preserve invalid acquisition rows in raw `Experiment` values and report them through validation issues | Accepted | Row removal, cropping, and repair are preprocessing decisions belonging to M3 |
| D023 | Encode dynamic runs as calibration data and all ramp-and-hold runs as held-out validation data in the dataset catalog | Accepted | Prevents accidental fitting on the predictive-validation programs |
| D024 | Emit deterministic Markdown and TOML audit artifacts from one command | Accepted | Supports both human review and machine regression checks without a tabular dependency |
| D025 | Estimate conversion-reference masses with robust local-linear fits at exact temperatures in physical windows | Accepted | Avoids dependence on whichever noisy acquisition sample happens to be nearest an endpoint |
| D026 | Use cubic local-polynomial conversion and heating-rate fits in a ±5 K window as `m3_default_v1` | Accepted | Real-data sensitivity reduces derivative roughness by about 47-fold versus raw finite differences while preserving peak height and reconstruction consistency |
| D027 | Diagnose conversion bounds, reversals, and negative derivatives without clipping or monotonic repair | Accepted | Keeps measurement and preprocessing limitations visible for M4 uncertainty analysis |
| D028 | Keep CairoMakie in the documentation environment rather than the numerical runtime | Accepted | Preserves pure numerical functions and avoids imposing plotting dependencies on library users |
| D029 | Carry 120/150 °C and 650/700/750/highest-common endpoint profiles into M4 sensitivity analysis | Accepted | Endpoint choice shifts conversion enough to matter but is not replicate uncertainty |
| D030 | Interpolate each experiment at the first acquisition-order upward crossing and prohibit extrapolation or silent monotonic repair | Accepted | Preserves measured ordering while exposing multiple crossings and unsupported conversion points |
| D031 | Canonically order experiments by identifier and fingerprint source checksums plus preprocessing and M4 configuration | Accepted | Makes numerical results and provenance invariant to caller order |
| D032 | Use ordinary least squares with transformed two-sided Student-t slope intervals for Friedman, KAS, FWO, and Starink | Accepted | Gives the linear-method intervals an explicit statistical meaning with four independent heating programs |
| D033 | Flag linear points below a configurable R² threshold without deleting their estimates | Accepted | Separates computational validity from weak straight-line support and keeps the 25% mixture tail visible |
| D034 | Minimize the true advanced-Vyazovkin pairwise objective using log-scaled trapezoidal integrals and a bounded scan plus golden-section refinement | Accepted | Removes the legacy hard-coded target of 12 and supports arbitrary experiment counts without numerical underflow |
| D035 | Use the Vyazovkin–Wight pairwise-variance Fisher procedure for advanced-method confidence intervals | Accepted | Separates statistical uncertainty from optimizer objective quality and reports bound truncation explicitly |
| D036 | Retain FWO as a compatibility method and prefer Starink among linear integral approximations | Accepted | Preserves thesis comparability while acknowledging the coarser Doyle approximation used by FWO |
| D037 | Report conversion-endpoint sensitivity separately from regression and Fisher intervals | Accepted | Changing reference temperatures changes the conversion coordinate and is not sampling uncertainty |
| D038 | Use one Fraser–Suzuki convention with a stable Gaussian limit and continuous zero outside the logarithm domain | Accepted | Removes singular evaluation near zero skew and makes invalid-domain behavior explicit |
| D039 | Fit measured `dalpha/dT` between exact first-upward 0.05 and 0.95 conversion crossings without clipping negative observations | Accepted | Keeps M5 aligned with preprocessing and preserves visible measurement limitations |
| D040 | Enforce peak constraints with smooth transformations and deterministic multistart L-BFGS | Accepted | Prevents label switching and order-dependent sequential bound tightening while exposing local minima |
| D041 | Retain raw BIC/AICc minima but require convergence, local identifiability, inactive bounds, and minimum component area for structural selection | Accepted | A lower residual criterion alone does not make an extra component numerically defensible |
| D042 | Report local physical-parameter Student-t intervals alongside Jacobian rank, normalized condition, correlation, and residual diagnostics | Accepted | Separates conditional uncertainty from optimizer quality while making non-identifiability visible |
| D043 | Keep shared-skew joint multi-rate fitting as a diagnostic, not the M5 default | Accepted | Two composition fits do not converge and two exceed the predeclared 5% objective-loss gate |
| D044 | Assign only ordered peak numbers in M5 and no chemical identities | Accepted | Fraser–Suzuki fit quality alone cannot identify a reaction or feedstock constituent |
| D045 | Keep M6 at the overall-conversion level and do not propagate M5 peaks | Accepted | M5 did not establish one cross-rate peak count or a joint component identity that passes its gate |
| D046 | Fit `ln(A/min^-1)` as the only intercept in a QR-based log-rate design | Accepted | Removes the duplicated intercept/amplitude non-identifiability and normal-equation instability in the legacy fits |
| D047 | Register 17 fixed equations plus separate two- and three-exponent Šesták–Berggren candidates | Accepted | Preserves the thesis model ideas while making every domain and parameter count explicit |
| D048 | Select with BIC after rank, condition, energy, and exponent screening, then require leave-one-heating-rate-out adequacy | Accepted | In-sample residual improvement alone is insufficient evidence of predictive kinetics |
| D049 | Report local Student-t coefficient intervals separately from correlation, autocorrelation, criterion ambiguity, and cross-validation | Accepted | The densely sampled derivative points do not justify treating conditional covariance as complete uncertainty |
| D050 | Treat the same-data compensation line as a diagnostic evaluated at an independent reference energy | Accepted | Strong `ln A`–`E` correlation can be mathematical or error-induced and does not establish a physical law |
| D051 | Reject a constant kinetic triplet as the default real-data predictive model | Accepted | All five selected composition models exceed the predeclared held-out-rate log-RMSE limit of 0.25 |

New decisions must be appended rather than silently changing prior entries. Superseded decisions should point to their replacements.

## 13. Risk register

| Risk | Impact | Mitigation |
|---|---|---|
| Some auxiliary ramp-and-hold channel labels are provisional | Incorrect DSC or gas-channel interpretation | Restrict validation to confirmed temperature/time/mass/segment columns until M2 resolves them |
| Exact nitrogen grade and historical correction provenance are unavailable | Incomplete reproducibility metadata | Report them as unknown and preserve all raw channels and provenance |
| Only a few heating rates are available | Wide or fragile uncertainty estimates | Report limitations; use method-specific intervals and sensitivity analyses |
| Peak parameters are non-identifiable | Unstable chemical interpretation | Ordered constraints, profile/bootstrapped uncertainty, model selection, joint fitting |
| Smoothing dominates derivative results | Biased Friedman and deconvolution estimates | Physical-window tuning and preprocessing sensitivity report |
| Legacy workspaces mix raw and derived variables | Accidental circular validation | Load only named raw variables from a documented allowlist |
| Optimizers converge to boundaries/local minima | Misleading parameter estimates | Bounded solvers, multistart diagnostics, convergence and boundary flags |
| Large binary artifacts overwhelm version control | Slow clones and unclear provenance | Keep generated outputs ignored; define a deliberate data-storage strategy |
| Model states leave `[0,1]` | Complex or nonphysical predictions | Bound-aware formulation, callbacks, and explicit solver failure |
| Compensation effect creates circular inference | Overconfident `A` estimates | Separate estimation stages and test identifiability/predictive value |
| M2 imports canonical MAT snapshots rather than reparsing NETZSCH text exports | Snapshot/export divergence could go unnoticed after source changes | Retain every raw-export path, verify MAT hashes, preserve the M0 exact-comparison evidence, and add direct text parsing if the exports change |
| MAT.jl brings a relatively large HDF5/binary transitive dependency graph | Slower first instantiation and a larger environment | Keep MAT.jl as the only M2 third-party runtime dependency and lock all resolved versions |
| Some real-data endpoint profiles and high-conversion linear fits are highly sensitive or weakly linear | Overinterpretation of apparent activation energy, especially for the 25% and 50% mixtures | Retain R² and optimizer warnings, report endpoint sensitivity separately, and avoid mechanistic claims from isolated profile excursions |
| Deconvolution residuals are serially correlated and many parameters are highly correlated | BIC/AICc and local covariance intervals can appear more decisive than warranted | Use one common deterministic criterion grid, retain Durbin–Watson/correlation diagnostics, and treat intervals as conditional |
| Selected real-data peak counts vary across composition and heating rate | Forced cross-rate labels could create fictitious component kinetics | Keep peaks ordered but chemically unnamed; require an explicit matching gate before M6 component propagation |
| Shared-skew joint fits fail some real-data convergence/objective gates | A convenient common-component model could degrade the data description | Retain joint fitting as a sensitivity diagnostic and use independent fits by default |
| Empirical M6 exponents approach configured plausibility bounds and have very high parameter correlations | Apparently precise reaction-model labels could be unstable or nonphysical | Retain bounds/correlation warnings, use labels only as empirical equations, and require prediction rather than mechanistic interpretation |
| Every selected constant-triplet model fails cross-rate prediction | M7 simulations based on one triplet could reproduce fit data yet fail new temperature programs | Keep the constant fits as frozen baselines and evaluate conversion-dependent kinetics against all held-out ramp-plus-hold runs |
| Compensation lines have high R² but require substantial extrapolation to some reference energies | Extrapolated pre-exponential factors could look authoritative despite weak support | Report the fitted energy range and interval warnings; never replace direct multi-rate estimates with compensation values |

## 14. Definition of done for every task

A checkbox is completed only when:

- implementation is committed to the intended module
- public behavior is documented
- relevant tests pass
- units and failure modes are explicit
- a real or synthetic example demonstrates the behavior
- generated evidence is linked from this plan or the appropriate validation document
- no unrelated legacy file was modified

## 15. Progress-update protocol

At the end of each work session:

1. Update the milestone table and affected checkboxes.
2. Add links to tests, reports, or figures that demonstrate completion.
3. Record new decisions in the decision log.
4. Add newly identified risks to the risk register.
5. Note blockers under the relevant milestone rather than marking incomplete work as complete.
6. Update the `Last updated` date.

## 16. Scientific references to anchor implementation

- Vyazovkin, S. et al. (2011), “ICTAC Kinetics Committee recommendations for performing kinetic computations on thermal analysis data,” *Thermochimica Acta* 520, 1–19. <https://doi.org/10.1016/j.tca.2011.03.034>
- Vyazovkin, S. (2001), “Modification of the integral isoconversional method to account for variation in the activation energy,” *Journal of Computational Chemistry* 22, 178–183. <https://doi.org/10.1002/1096-987X(20010130)22:2%3C178::AID-JCC5%3E3.0.CO;2-%23>
- Vyazovkin, S. and Wight, C. A. (2000), “Estimating realistic confidence intervals for the activation energy determined from thermoanalytical measurements,” *Analytical Chemistry* 72, 3171–3175. <https://doi.org/10.1021/ac000210u>
- Starink, M. J. (2003), “The determination of activation energy from linear heating rate experiments: a comparison of the accuracy of isoconversion methods,” *Thermochimica Acta* 404, 163–176. <https://doi.org/10.1016/S0040-6031(03)00144-8>
- Fraser, R. D. B. and Suzuki, E. (1969), “Resolution of overlapping bands: Functions for simulating band shapes,” *Analytical Chemistry* 41, 37–39. <https://doi.org/10.1021/ac60270a007>
- Šesták, J. and Berggren, G. (1971), “Study of the kinetics of the mechanism of solid-state reactions at increasing temperatures,” *Thermochimica Acta* 3, 1–12. <https://doi.org/10.1016/0040-6031(71)85051-7>
- Roduit, B. (2000), “Computational aspects of kinetic analysis. Part E: The ICTAC Kinetics Project—numerical techniques and kinetics of solid state processes,” *Thermochimica Acta* 355, 171–180. <https://doi.org/10.1016/S0040-6031(00)00447-0>
- Vyazovkin, S. (2021), “A time to search: finding the meaning of the kinetic compensation effect,” *Molecules* 26, 3077. <https://doi.org/10.3390/molecules26113077>
- Barrie, P. J. (2012), “The mathematical origins of the kinetic compensation effect: 1. the effect of random experimental errors,” *Physical Chemistry Chemical Physics* 14, 318–326. <https://doi.org/10.1039/C1CP22666E>

Additional method references and exact formula conventions are recorded in `docs/methods.md`.
