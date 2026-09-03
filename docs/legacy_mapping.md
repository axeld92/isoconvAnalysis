# Legacy Workflow and Artifact Map

Last updated: 2026-08-17  
M0 status: Complete; exploratory files remain archive-only

## 1. End-to-end workflow recovered from the MATLAB archive

```text
NETZSCH exports / MAT workspaces
        |
        v
crop, normalize mass, smooth, differentiate
        |
        +--> TGA/DTG diagnostics
        |
        +--> isoconversional E_alpha
        |       Friedman / KAS / FWO / Starink / advanced Vyazovkin
        |
        +--> Fraser-Suzuki DTG deconvolution
                |
                +--> component-specific E_alpha
                +--> compensation / pre-exponential estimates
                +--> Sestak-Berggren reaction-model fits
                +--> dynamic and isothermal ODE prediction

pure-component curves + mixture curve
        |
        v
additive reference and synergy difference
```

This workflow is the intellectual content to preserve. The two feedstocks are waste tire and *Acrocomia aculeata* endocarp, mixed on a dry-mass basis. Script ordering, global workspaces, hard-coded thresholds, and saved fit parameters are not part of the specification.

## 2. Dataset roles

| Artifact | Classification | Role in rewrite |
|---|---|---|
| `code/2021corridas/ExpDat_*wt*hr*.txt` | Core raw evidence | Primary dynamic source |
| `code/2021corridas/2021corridas.mat` | Core raw snapshot | Canonical MAT import/regression fixture source |
| `code/2021mar.mat` | Core alias snapshot | Exact `St01`–`St20` mapping evidence |
| `code/isotherm.mat` | Core, schema provisional | Isothermal validation candidate |
| `code/2021marcorrected.mat` | Reference | Legacy preprocessing comparison |
| `code/corridasjulio.mat`, `code/julio2.mat` | Reference | Final-analysis workspace inputs; not raw |
| `code/a*neum.mat`, `code/*_nuevo.mat` | Reference results | Legacy result comparison only |
| `code/analisisarticuloneumatics.mat`, `code/articulo/articuloneumaticos.mat` | Reference results | Publication workspace snapshots |
| `code/matlab.mat`, `code/oct2020*.mat`, `code/2020oct*.mat` | Superseded/reference | Earlier data generations |
| `code/TGAnalysis-master/`, `code/github_repo/` | Vendored/unrelated | Exclude from the Julia port |

## 3. Script classification

### 3.1 Core scientific ideas

These files contain concepts or equations that should be represented in the new implementation, after correction and validation.

| Area | Files |
|---|---|
| Data conditioning | `cleandata.m`, `cleandata2.m`, `cleandatam.m`, `cleandatami.m`, `cleandatami2.m`, `diff_alpha.m`, `removeerror.m`, `removeerror2.m`, `removeerrorm.m`, `getind.m`, `Tfit.m`, `timefit.m`, `smoothdadT.m`, `smoothdadT2.m` |
| Isoconversional analysis | `friedman.m`, `FWO.m`, `KAS.m`, `Starink.m`, `Starink3.m`, `Starinkv.m`, `Vyazovkin.m`, `Vyazovkina.m`, `Vyazovkinfor3.m`, `Vyazovkinv.m`, `friedmannea.m`, `phiofe2.m`, `phiofe3.m`, `innerint2.m`, `Jintegs.m`, `polyparci.m` |
| Peak model and deconvolution | `fs_function.m`, `fs_mixture.m`, `fs_mixture2.m`, `fs_fit.m`, `fs_out.m`, `fs_out2.m`, `fslogn.m`, `deconvolution.m`, `deconvolution4.m`, `deconvasym.m`, `deconvolve2.m`, `deconvolve4.m`, `deconvolveall.m`, `deconvolveg.m`, `peaks3or4.m` |
| Kinetic triplet / compensation | `comeffect.m`, `comeffect2.m`, `comeffect3.m`, `compeffect.m`, `compeffects.m`, `icompeffect.m`, `preexpcomp.m`, `preexpfactor.m`, `sb.m`, `sblin.m`, `sbrlin.m`, `kineticanalysis.m`, `lambdamasterplot.m` |
| Simulation and mixture interaction | `simulacionvyazov.m`, `simuisotherm.m`, `simsigmoid.m`, `simsigmoid2.m`, `simsigmoid3.m`, `movietga.m`, `synergy.m` |
| TGA characteristics | `tgaparam.m` |

Core means “preserve the idea,” not “translate the code line by line.”

### 3.2 Workflow and reporting references

These scripts show which combinations were run or how results were presented:

`calc_EA.m`, `calc_EA2.m`, `isoconvanalysis.m`, `deconvanalysis.m`, `deconvolve.m`, `ea_deconv.m`, `ea_deconv4.m`, `ea_deconvv.m`, `ea_deconvv4.m`, `integrate_all_dadt.m`, `integrate_all_dadTT.m`, `integratedadt.m`, `kinectideconvolutionanalysis.m`, `modelfitting.m`, `exportcleandata.m`, `getdadT.m`, `plot_ind_comp.m`, `plotafa.m`, `plotalpha.m`, `plotcompenseffect.m`, `plotconversion.m`, `plotdadt.m`, `plotfiguresforthesis.m`, `plotfsfunctions.m`, `plotmanddm.m`, `plotmassloss.m`, `plotmasslossi.m`, `plottingforthesis.m`, `showfsdeconv.m`, and `scriptgenera.m`.

They are useful for reconstructing expected outputs but are not acceptable entry points for the Julia package.

### 3.3 Exploratory, duplicated, or abandoned work

These files may contain useful clues but are not specification sources:

`Copy_of_ea_deconv.m`, `Copy_of_prueba_fs.m`, `Untitled.m`, `aicvy2.m`, `cargardatosoct2020.m`, `cargardatossept2020.m`, `comparison.m`, `eavyazovprueba.m`, `forreddit.m`, `forredditiguesss.m`, `friedman2.m`, `lookingforproblems.m`, `newapproach.m`, `probandoprobando.m`, `probarlinreg.m`, `prueba_fs.m`, `pruebaaaaaa.m`, `pruebadatospubli.m`, `pruebaglobalsearch.m`, `seguro que algo de dmd.m`, `somethingsomething.m`, `sortfilas.m`, `testnewapproach.m`, `toofasttoofourier.m`, `toyprob.m`, `toyprob2.m`, `toyprob3.m`, `trysmoothdata.m`, and `varofe.m`.

`stylizeddoconvolution.txt` is a later prose/code sketch and is treated as a design reference, not historical ground truth.

### 3.4 Unrelated demonstrations or third-party material

- `CK_guillermo.m`, `Comp_Koopman.m`, `collatz.m`, `datadrivenscienceandengineering.m`, and `svdanalysis.m` are unrelated numerical demonstrations.
- `code/TGAnalysis-master/` and `code/github_repo/` are duplicate copies of a third-party TGAnalysis project.
- Archives such as `github_repo.zip`, `TGAnalysis-master.zip`, and `resultados.rar` are not inputs to the rewrite.

## 4. Expected scientific outputs recovered from artifacts

| Output family | Legacy evidence | Rewrite status |
|---|---|---|
| TGA and DTG curves by composition and heating rate | `resultados/tgadtg/`, `articulo/tgadtg.*`, `plotmassloss.m` | Required |
| Conversion curves | `resultados/marzo/`, `resultados/octubre/`, `plotconversion.m` | Required QC output |
| Activation energy versus conversion for five compositions | `resultados/analisis isoconv/`, `articulo/EAs.pdf`, `calc_EA2.m` | Required |
| Fraser–Suzuki fits for each composition/rate | `resultados/deconv/`, `articulo/{5K,10K,15K,20K} FS.*` | Required after model-selection validation |
| Component-specific kinetic parameters | `ea_deconv*.m`, `resultados/A/` | Required after deconvolution validation |
| Compensation-effect plots | `articulo/compeffect.*`, `plotcompenseffect.m` | Candidate; must pass identifiability review |
| Sesták–Berggren fits and intervals | `articulo/SB*`, `sbrlin.m` | Candidate; must pass predictive validation |
| Dynamic simulations | `resultados/simu*`, `articulo/simus.*` | Required predictive output |
| Ramp-and-hold predictions at 450/500/550 °C | `resultados/isot/`, `simuisotherm.m`, `code/isotherm.mat` | Required held-out validation; no parameter refitting |
| Synergy curves for 25/50/75 wt% mixtures | `resultados/synergy/`, `synergy.m` | Required descriptive output |

No legacy plot is a pixel-matching target. Release acceptance is based on the scientific content of the outputs, especially experimental-versus-predicted conversion, the imposed temperature program, ramp/hold residuals, and errors by composition and hold temperature.

## 5. Known legacy incompatibilities that affect mapping

- Raw direct-import arrays use mass in column 4, while later five-column arrays use mass in column 3.
- Four-column processed arrays use column 3 for either `dalpha/dT`, `dalpha/dt`, or a raw signal depending on the producing function.
- Some routines add 273.15 to already-Kelvin processed temperature.
- `cleandata.m` hard-codes 150–700 °C and filters derivatives with unexplained numerical cutoff values.
- Several isoconversional functions accept only four runs even though the underlying method is general.
- `deconvolve.m` ignores its input curve and fits the hard-coded `St05` workspace variable.
- The legacy four-peak area-fraction calculation omits peak four and repeats peak three.
- `deconvolveall.m` contains undefined and misassigned variables, so it is not an executable specification.
- Deconvolution workspaces contain both independent fits and sequentially tightened-bound fits,
  so saved results can depend on experiment order.
- `kineticanalysis.m` places two constant columns in the same design matrix, so its
  pre-exponential intercept and extra amplitude coefficient are not separately identifiable;
  it also solves least squares through normal equations.
- `sblin.m` repeats the same amplitude/pre-exponential confounding after fixing `A`, while
  `sbrlin.m` removes the intercept but reports observation-sized prediction widths as though
  they were parameter confidence intervals.
- `sb.m` searches arbitrary `[-100,100]` empirical exponents without a predictive gate, and
  `preexpcomp.m` is syntactically incomplete.
- `kinectideconvolutionanalysis.m` depends on globals, hard-codes `St16`, and leaves component
  states unconstrained; its component propagation is not an executable specification.
- `compeffects.m` infers compensation from a single heating program and contains an incorrect
  Ginstling–Brounshtein expression. A high same-data compensation correlation is retained only
  as a diagnostic in M6.

The Julia import layer will therefore identify schemas by named source/configuration, never by column count alone.

M6 replaces these kinetic paths with one named registry, a single-intercept QR design,
conditional coefficient intervals, structural model screening, and leave-one-heating-rate-out
validation. The real-data gate rejects every selected constant-triplet candidate as a default
predictive model; no M5 peak is propagated into a reaction identity.

## 6. M0 conclusions

- The main dynamic design and `St01`–`St20` mapping are recovered with exact evidence.
- The authoritative raw column meanings and units are recovered from instrument headers.
- The scientific workflow and expected output families are recovered.
- Sample identities, dry-mass composition basis, nitrogen atmosphere, 150–700 °C primary conversion convention, and ramp-and-hold validation role are confirmed.
- Exact nitrogen grade and some historical export/correction provenance remain non-blocking metadata gaps.
