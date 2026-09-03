# Contributing

Use Julia 1.11 and keep numerical work separate from plotting and file I/O. Public functions
need docstrings, new behavior needs tests, and scientific assumptions must be reflected in
the masterplan or method documentation.

Before submitting a change, run:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project=docs scripts/format.jl --check
julia --startup-file=no --project=docs docs/make.jl
```

To apply the formatter, omit `--check`. Do not modify files under `code/`; they are legacy
evidence for comparison and provenance.
