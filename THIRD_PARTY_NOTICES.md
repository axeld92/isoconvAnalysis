# Third-party notices

Last reviewed: 2026-09-05

## TGAnalysis archive

The legacy evidence directory retains data, results, and attribution from two former
copies of the external TGAnalysis project:

- `code/github_repo/`
- `code/TGAnalysis-master/TGAnalysis-master/`

Upstream project: [mayankskii/TGAnalysis](https://github.com/mayankskii/TGAnalysis)

TGAnalysis is Copyright (c) 2019 Mayank Mahawar and distributed under the MIT License. Each
archived copy retains the upstream `LICENSE` and `README.md`. The two local license files are
byte-identical and have SHA-256:

```text
9f64ce19c811f17124ac3bc679b069d51748eab537f2290a33171e2d76688966
```

The TGAnalysis MATLAB source was removed from the project tree. Preserve the copyright
notice and MIT permission notice if the retained material is redistributed.

## Julia dependencies

Julia packages are resolved by `Pkg` from the locked manifests and are not vendored into this
repository. Direct-dependency license evidence is recorded in
[`docs/dependencies.md`](docs/dependencies.md). Binary artifacts downloaded by Julia packages
may carry their own upstream notices; those notices must be collected separately if a future
release redistributes a compiled environment rather than source, manifests, and generated
figures only.
