# Development

From the repository root, instantiate and test the runtime environment:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

Instantiate the isolated documentation/tooling environment, check formatting, and build the
manual:

```sh
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=docs scripts/format.jl --check
julia --startup-file=no --project=docs docs/make.jl
```

Regenerate the M2 raw-data audit, M3 preprocessing evidence, M4 isoconversional evidence,
M5 deconvolution evidence, and M6 kinetic-triplet evidence with:

```sh
julia --startup-file=no --project=. scripts/audit_legacy_data.jl
julia --startup-file=no --project=docs scripts/audit_preprocessing.jl
julia --startup-file=no --project=docs scripts/audit_isoconversional.jl
julia --startup-file=no --project=docs scripts/audit_deconvolution.jl
julia --startup-file=no --project=docs scripts/audit_reaction_models.jl
```

Run `scripts/format.jl` without `--check` to apply formatting. Dependency rationale and
license information are recorded in the repository file `docs/dependencies.md`.
