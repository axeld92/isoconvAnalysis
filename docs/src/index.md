# IsoconversionalAnalysis.jl

`IsoconversionalAnalysis.jl` is a reproducible Julia rewrite of a master-thesis
thermogravimetric and isoconversional-analysis workflow.

The project is being built milestone by milestone. M0–M6 now provide the scientific
specification, reproducible Julia foundation, checksum-verified data ingestion, validated
preprocessing and isoconversional-analysis pipelines, and constrained Fraser–Suzuki
deconvolution, plus an audited overall-conversion kinetic-triplet and reaction-model layer.
Scientific requirements, terminology, and decisions are tracked in the
repository file `MASTERPLAN.md`.

```@docs
IsoconversionalAnalysis
```

## Foundation API

```@autodocs
Modules = [IsoconversionalAnalysis]
Order = [:type, :function]
```
