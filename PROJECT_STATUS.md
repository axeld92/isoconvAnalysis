# Project status and next steps

- Last updated: 2026-09-04
- Current milestone: M6 complete; M7 ready to start
- Repository: <https://github.com/axeld92/isoconvAnalysis>

## Project summary

`IsoconversionalAnalysis.jl` is a from-scratch Julia 1.11 implementation of the useful
scientific ideas from the original MATLAB master-thesis workflow. It analyzes
thermogravimetric experiments for waste tire, *Acrocomia aculeata* endocarp, and their dry-mass
mixtures. The rewrite prioritizes traceability, explicit units, typed results, visible failure
modes, synthetic validation, and prediction on experiments excluded from parameter fitting.

The original `code/` directory remains an immutable historical archive. Legacy results are
comparison evidence, not numerical ground truth, and known implementation defects are not
reproduced in the Julia package.

## Experimental and scientific contract

- Mixture percentages are dry-mass fractions. Samples were dried at 105 °C for at least 48 h.
- The confirmed atmosphere is nitrogen in the purge and protective channels; the exact
  historical gas grade remains unknown.
- The primary conversion coordinate uses masses evaluated at 150 °C and 700 °C. Endpoint
  sensitivity profiles remain separate from statistical confidence intervals.
- Twenty dynamic runs cover five compositions at nominal heating rates of 5, 10, 15, and
  20 K/min and form the calibration dataset.
- Fifteen ramp-and-hold runs in `code/isotherm.mat` are reserved for out-of-sample validation.
  They must never be used for parameter fitting or model selection.
- Reaction-model and deconvolution labels are empirical descriptions unless independent
  evidence establishes a physical or chemical identity.

## Completed work

| Milestone | Delivered capability | Principal evidence |
|---|---|---|
| M0 | Scientific requirements, data dictionary, units, experiment map, and legacy classification | [`docs/data_dictionary.md`](docs/data_dictionary.md), [`docs/methods.md`](docs/methods.md) |
| M1 | Reproducible Julia package, locked environments, CI, formatting, tests, and documentation | [`docs/m1_foundation.md`](docs/m1_foundation.md) |
| M2 | Checksum-verified MAT ingestion for all 35 configured experiments and deterministic data audits | [`docs/data_audit.md`](docs/data_audit.md) |
| M3 | Source-traceable preprocessing, conversion calculation, smoothing/derivative sensitivity, and QC figures | [`docs/preprocessing_audit.md`](docs/preprocessing_audit.md) |
| M4 | Friedman, KAS, FWO, Starink, and advanced Vyazovkin methods with documented uncertainty | [`docs/isoconversional_audit.md`](docs/isoconversional_audit.md) |
| M5 | Constrained Fraser–Suzuki fits, peak-count comparison, identifiability diagnostics, and joint-fit evaluation | [`docs/deconvolution_audit.md`](docs/deconvolution_audit.md) |
| M6 | Nineteen reaction models, identifiable kinetic-triplet regression, compensation diagnostics, and cross-rate validation | [`docs/reaction_model_audit.md`](docs/reaction_model_audit.md) |

The complete package suite currently passes 369 tests, including Aqua package-quality checks.
JuliaFormatter and Documenter gates pass, and the generated M2–M6 audit artifacts are
deterministic across repeated runs.

## Important findings and accepted decisions

1. The validated M3 preprocessing profile uses cubic local-polynomial fits in a ±5 K window.
   Negative derivatives and conversion reversals remain visible rather than being clipped or
   silently repaired.
2. The five M4 isoconversional methods operate on a shared first-upward-crossing conversion
   grid. Advanced Vyazovkin minimizes the actual pairwise objective and is not hard-coded for
   four experiments.
3. Independent M5 peak fits remain the default. Real experiments select different peak counts,
   and the evaluated joint shared-skew model does not pass every convergence and objective-loss
   gate. Peaks therefore have ordered numbers but no forced cross-rate or chemical identity.
4. M6 removes the legacy duplicated-intercept problem by fitting exactly one `ln(A/min^-1)`
   intercept. It retains all candidate results, structural diagnostics, local conditional
   intervals, and leave-one-heating-rate-out errors.
5. Every selected constant-triplet model fails the predeclared M6 cross-rate log-RMSE limit of
   0.25. Constant triplets are consequently descriptive benchmarks, not accepted predictive
   mechanisms.
6. Compensation-effect fits have high apparent correlations but require extrapolation in all
   20 runs. They remain diagnostics and are never used to replace direct multi-rate estimates.

## Repository layout

```text
src/       Julia numerical implementation
test/      synthetic, failure-mode, and real-source integration tests
config/    dataset catalog and frozen analysis defaults
scripts/   deterministic audit and development commands
docs/      methods, reports, machine-readable evidence, figures, and manual sources
code/      immutable legacy MATLAB scripts, data, and historical outputs
```

The detailed scientific plan, decision log, risks, and milestone acceptance criteria remain in
[`MASTERPLAN.md`](MASTERPLAN.md). This file is the concise operational snapshot.

## Immediate next milestone: M7 simulation and predictive validation

M7 must determine whether the validated conversion-dependent isoconversional information can
predict temperature programs that differ from the calibration ramps. Work should proceed in
the following order:

1. Define typed linear-ramp, isothermal-hold, and ramp-plus-hold temperature programs with
   explicit kelvin/minute units and validation.
2. Define the simulation configuration and typed result/diagnostic contracts before selecting
   an ODE solver.
3. Implement the M6 constant-triplet formulation as a frozen benchmark without refitting it.
4. Construct an internally consistent conversion-dependent `E(alpha)` and `A(alpha)` path from
   calibration data only. Define interpolation support and boundary behavior explicitly.
5. Enforce finite physical states with `0 <= alpha <= 1`; terminate with a diagnostic instead
   of hiding invalid or complex states.
6. Validate linear-ramp and ramp-plus-hold integration on synthetic cases with known solutions,
   irregular sampling, and solver-failure fixtures.
7. Predict all 15 ramp-and-hold experiments without refitting. Report conversion overlays,
   measured/model temperature programs, and residuals separately for ramp and hold phases.
8. Compare the constant-triplet and conversion-dependent approaches using frozen error and
   physical-validity gates. Retain negative results explicitly.

### M7 exit evidence

M7 is complete only when:

- every prediction is reproducible, finite, and accompanied by solver diagnostics;
- the physical conversion domain is enforced by design;
- no held-out ramp-and-hold observation influences parameter estimation or model choice;
- residuals and summary errors are reported by composition, hold temperature, ramp, and hold;
- synthetic recovery and all existing package tests pass; and
- the preferred predictive formulation is selected from held-out evidence or all candidates
  are explicitly rejected.

## Later roadmap

### M8 — Synergy and reporting

- Implement dry-mass additive references for mass, conversion, and conversion rate.
- Quantify mixture deviations with uncertainty bands and clearly defined signs.
- Produce the standard scientific figures and one reproducible reference report.
- Regenerate the analysis and provenance record from one documented command.

### M9 — Release

- Complete the corrected Julia-versus-legacy comparison and document intentional deviations.
- Finish user/developer documentation and review data-sharing and dependency licenses.
- Run the clean-environment validation matrix and inspect CI on GitHub.
- Select a project license, tag version 1.0, and archive the validation bundle.

## Reproduction commands

From the repository root:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project=docs scripts/format.jl --check
julia --startup-file=no --project=docs docs/make.jl
```

The individual M2–M6 audit commands are listed in
[`docs/src/development.md`](docs/src/development.md).
