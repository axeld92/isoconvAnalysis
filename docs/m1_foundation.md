# M1 Julia foundation report

Date: 2026-08-17  
Status: complete locally

## Outcome

M1 establishes a loadable Julia 1.11 package with no third-party runtime dependencies.
Configuration, logging, testing, formatting, documentation, and CI are present before any
scientific numerical implementation begins.

The environment is split by responsibility:

| Environment | Files | Direct dependencies | Exact verified versions |
|---|---|---|---|
| Runtime | `Project.toml`, `Manifest.toml` | Julia `Logging` and `TOML` standard libraries | Julia 1.11.6, Logging 1.11.0, TOML 1.0.3 |
| Test | `test/Project.toml`, `test/Manifest.toml` | package under test, `Test`, `Logging`, Aqua.jl | Aqua 0.8.16 |
| Documentation/tooling | `docs/Project.toml`, `docs/Manifest.toml` | package, Documenter.jl, JuliaFormatter.jl | Documenter 1.17.0, JuliaFormatter 2.12.5 |

This separation prevents developer tools from becoming runtime dependencies while keeping
all three environments reproducible.

## Implemented foundation

- `IsoconversionalAnalysis` package entry point and public foundation API.
- Typed `AnalysisConfig`, `ProjectConfig`, `UnitConfig`, and `ConversionConfig` values.
- TOML parsing with explicit errors for missing fields, invalid types, non-finite endpoints,
  unsupported log levels, noncanonical unit labels, and invalid conversion intervals.
- Configuration-relative path resolution, independent of the caller's working directory.
- Standard Julia structured logging with a configurable minimum level.
- Default M0 decisions in `config/analysis_defaults.toml`, including the primary 150–700 °C
  conversion interval and endpoint sensitivity cases.
- Unit tests and Aqua package-quality checks.
- Blue-style JuliaFormatter command and strict Documenter build.
- GitHub Actions jobs for Julia 1.11 tests, formatting, and documentation.
- Ignore rules for generated results, derived data, built documentation, local archives, and
  large legacy video artifacts. Legacy source and scientific MAT inputs remain visible.

## Verification evidence

The following commands passed from the repository root on Julia 1.11.6:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. -e 'using IsoconversionalAnalysis; load_config("config/analysis_defaults.toml")'
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project=docs scripts/format.jl --check
julia --startup-file=no --project=docs docs/make.jl
```

The test suite passes 31 checks. Aqua covers method ambiguities, unbound type parameters,
undefined exports, stale dependencies, compatibility bounds, method piracy, and persistent
tasks. The documentation build includes doctests, cross-reference validation, and exported
docstring coverage.

The CI workflow contains these same gates. It has not run on a remote CI service because the
local repository currently has no configured Git remote; this does not affect local
instantiation or the fresh-clone command sequence.

## Deferred deliberately

- Data structures and MAT-file dependencies begin in M2.
- A physical-units package will be evaluated when numerical arrays are introduced; M1 uses
  typed field names, canonical labels, and boundary validation.
- Machine-readable log sinks can be added when batch pipelines need them; scientific code
  already emits backend-independent structured metadata.
- The implementation's own distribution license remains a user decision before public
  release. Dependency licenses are recorded in `docs/dependencies.md`.
