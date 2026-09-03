# Dependency policy and licenses

Last reviewed: 2026-08-17

The package runtime is kept deliberately small. Numerical and data-format packages are added
only when their milestone supplies a concrete requirement, validation plan, and maintenance
rationale. M2 adds MAT.jl as the first third-party runtime dependency; M4 adds Distributions.jl
for method-specific confidence quantiles; M5 adds Optim.jl for deterministic constrained
nonlinear fitting. M6 reuses the standard-library linear algebra and Distributions.jl
facilities and adds no runtime dependency. Compatible versions
are constrained in the relevant `Project.toml`; the locked manifests record exact resolved
versions.

| Dependency | Scope | Choice and rationale | License |
|---|---|---|---|
| Julia 1.11 | Runtime | Current project baseline with stable package environments and modern language features. | MIT |
| `Logging` | Runtime, standard library | Structured event metadata without a third-party logging dependency. | MIT (Julia distribution) |
| `TOML` | Runtime, standard library | Native parsing for human-readable, version-controlled configuration. | MIT (Julia distribution) |
| `Dates` | Runtime, standard library | Typed acquisition timestamps from catalog metadata. | MIT (Julia distribution) |
| `SHA` | Runtime, standard library | Source-file fingerprints checked before MAT variables are loaded. | MIT (Julia distribution) |
| `Statistics` | Runtime, standard library | Median sampling intervals and OLS audit diagnostics. | MIT (Julia distribution) |
| `Printf` | Runtime, standard library | Deterministic numeric formatting in the Markdown audit. | MIT (Julia distribution) |
| `LinearAlgebra` | Runtime, standard library | Stable least-squares calculations for robust references, local polynomials, and fitted heating rates. | MIT (Julia distribution) |
| `MAT.jl` 0.12.1 | Runtime | Explicit named-variable access to the legacy MATLAB v5 sources; avoids importing complete workspaces. | MIT |
| `Distributions.jl` 0.25.130 | Runtime | Student-t and F quantiles for documented linear-regression and Vyazovkin–Wight confidence intervals. | MIT |
| `Optim.jl` 2.2.1 | Runtime | L-BFGS with finite-difference gradients for deterministic multistart Fraser–Suzuki fitting; physical constraints are encoded by smooth transformations. | MIT |
| `Aqua.jl` 0.8.16 | Test only | Detects package hygiene problems, dependency issues, ambiguities, and method piracy. | MIT |
| `JuliaFormatter.jl` 2.12.5 | Development/docs environment | Deterministic Blue-style formatting and a CI-checkable style gate. | MIT |
| `Documenter.jl` 1.17.0 | Development/docs environment | Standard Julia API/manual generation with docstring coverage checks. | MIT |
| `CairoMakie.jl` 0.15.13 | Development/docs environment | Publication-quality static M3 QC figures without adding plotting to numerical package code. | MIT |

GitHub Actions used by CI are build infrastructure, not package dependencies. The workflow
pins their major release lines so routine upstream fixes remain available.

The new implementation itself does not yet declare a license. A project license must be
selected explicitly before public distribution; dependency licenses do not determine it.
