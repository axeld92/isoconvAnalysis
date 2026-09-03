# Scientific Method Contracts

Last updated: 2026-08-17  
M6 status: Kinetic-triplet and reaction-model contracts implemented and audited

## 1. Purpose

This is the equation registry for the Julia rewrite. It records the scientific ideas to preserve and fixes the notation and units that were inconsistent in the legacy MATLAB code.

These are implementation contracts, not claims that every method is appropriate for every dataset. Synthetic and predictive validation remain mandatory.

## 2. Common notation and units

| Symbol | Meaning | Canonical unit |
|---|---|---|
| `alpha` | conversion extent | dimensionless, normally 0–1 |
| `t` | elapsed time | min at the data boundary; converted consistently if a solver uses seconds |
| `T` | absolute temperature | K |
| `beta = dT/dt` | heating rate | K/min |
| `E_alpha` | apparent activation energy at conversion `alpha` | kJ/mol in public results |
| `A_alpha` | apparent pre-exponential factor | `min^-1` for the implemented dimensionless differential models |
| `R` | molar gas constant | `0.00831446261815324 kJ mol⁻¹ K⁻¹` when `E` is kJ/mol |
| `f(alpha)` | differential reaction model | dimensionless under the selected rate convention |

The base rate equation is

```text
dalpha/dt = A_alpha * exp(-E_alpha / (R*T)) * f(alpha)
```

For a linear ramp:

```text
dalpha/dt = beta * dalpha/dT
```

All logarithms below are natural logarithms unless explicitly written as `log10`.

## 3. Conversion definition

The legacy work used

```text
alpha(T) = (m_initial - m(T)) / (m_initial - m_final)
```

Mixtures were prepared on a dry-mass basis after oven drying at 105 °C for at least 48 hours. The accepted primary analysis coordinate is

```text
alpha_700(T) = [m(150 °C) - m(T)] / [m(150 °C) - m(700 °C)]
```

The two reference masses must be evaluated at the exact reference temperatures with a validated local robust fit or shape-preserving interpolation, not selected as arbitrary neighboring rows. The primary results must be accompanied by sensitivity calculations using a 120 °C initial reference and final references of 650 °C, 750 °C, and the highest common observed temperature without extrapolation.

This convention defines the conversion coordinate for kinetic comparison. It does not assert that decomposition is chemically complete at 700 °C.

## 4. Isoconversional methods

At fixed `alpha`, the reaction model term is constant across heating programs. Every method therefore evaluates different experiments at the same conversion, never merely at the same temperature or time.

### 4.1 Friedman differential method

```text
ln[(dalpha/dt)_(alpha,i)]
    = ln[A_alpha * f(alpha)] - E_alpha/(R*T_(alpha,i))
```

Regress `ln(dalpha/dt)` against `1/T`. The slope is `-E_alpha/R`.

For linear heating, `dalpha/dt` may be calculated as `beta * dalpha/dT`, provided both terms use compatible time and temperature units.

Legacy mapping: `friedman.m`. The implementation there contains an interpolation defect in the fourth derivative input and a parenthesis defect in heating-rate estimation; neither is part of the method contract.

### 4.2 Kissinger–Akahira–Sunose (KAS)

```text
ln[beta_i / T_(alpha,i)^2]
    = C_alpha - E_alpha/(R*T_(alpha,i))
```

Regress `ln(beta/T^2)` against `1/T`. The slope is `-E_alpha/R`.

Legacy mapping: `KAS.m`.

### 4.3 Flynn–Wall–Ozawa (FWO)

Natural-log convention used by the legacy code:

```text
ln(beta_i) = C_alpha - 1.052 * E_alpha/(R*T_(alpha,i))
```

Equivalent base-10 forms use a different coefficient. The code must bind the coefficient to the logarithm convention and test the pair together.

Legacy mapping: `FWO.m`.

### 4.4 Starink

```text
ln[beta_i / T_(alpha,i)^1.92]
    = C_alpha - 1.0008 * E_alpha/(R*T_(alpha,i))
```

Regress `ln(beta/T^1.92)` against `1/T`. The slope is `-1.0008 E_alpha/R`.

Legacy mapping: `Starink.m`, `Starink3.m`, and `Starinkv.m`.

### 4.5 Advanced Vyazovkin

For experiment `i` and a small conversion interval ending at `alpha`:

```text
J_i(E_alpha, alpha)
    = integral from t_(alpha-delta_alpha,i) to t_(alpha,i)
      of exp[-E_alpha/(R*T_i(t))] dt
```

Define the pairwise objective for `n` experiments:

```text
Phi(E_alpha) = sum_i sum_(j != i) J_i(E_alpha,alpha) / J_j(E_alpha,alpha)
```

Then

```text
E_alpha = argmin_E Phi(E)
```

When every `J_i` is identical, the objective equals `n*(n-1)`. This observation is a diagnostic; the algorithm must minimize `Phi` directly and must not minimize `abs(Phi - 12)`. The legacy value 12 only corresponds to four experiments and fails for other experiment counts.

Required configuration includes `delta_alpha`, energy bounds, integration rule/resolution, interpolation policy, and handling of objective boundary solutions.

Legacy mapping: `Vyazovkin.m`, `Vyazovkinv.m`, `Vyazovkinfor3.m`, `phiofe2.m`, `phiofe3.m`, `innerint2.m`, and `Jintegs.m`.

## 5. Uncertainty contracts

- Linear-method intervals are two-sided Student-t intervals for the OLS slope and its
  transformation to activation energy. R² is reported separately as a fit diagnostic.
- The advanced-Vyazovkin interval follows Vyazovkin and Wight. For `n` experiments, define

```text
S²(E) = [1/(n*(n-1))] * sum_i sum_(j != i) [J_i(E)/J_j(E) - 1]^2
Psi(E) = S²(E) / min_E S²(E)
```

  The confidence set contains energies for which
  `Psi(E) < F(confidence; n-1, n-1)`. Numerical roots on either side of the variance
  minimum provide the interval limits; a missing root is explicitly reported as truncation
  at the configured energy bound.
- Optimizer objective values and residuals are diagnostics, not confidence intervals.
- Preprocessing sensitivity and replicate variability must be reported separately from regression uncertainty.
- A result type must state confidence level, interval construction, assumptions, and number of independent experiments.

Legacy mappings include `polyparci.m` and the experimental interval calculation in `Vyazovkinv.m`. Those calculations are comparison targets only.

## 6. Fraser–Suzuki peak model

The legacy four-parameter peak is

```text
F(T; h, s, p, w)
  = h * exp{-[ln(2)/s^2] * [ln(1 + 2*s*(T-p)/w)]^2}
```

on the domain

```text
1 + 2*s*(T-p)/w > 0
```

with Gaussian limit

```text
F(T; h, 0, p, w)
  = h * exp{-4*ln(2)*[(T-p)/w]^2}
```

where `h` is peak height, `s` skew, `p` peak position, and `w` width. A mixture is a sum of peak functions.

Implementation requirements:

- stable evaluation around `s = 0`
- positive height and width
- explicit domain behavior
- ordered centers or another label-identifiability constraint
- derivative-rate and integrated-area reconstruction checks
- model-selection and uncertainty diagnostics

M5 freezes the following implementation contract:

- fit measured `dalpha/dT` on the exact first-upward crossings of `alpha = 0.05` and
  `alpha = 0.95`; do not clip negative observations
- constrain nonnegative height, bounded skew and width, and ordered centers separated by at
  least 8 K through smooth parameter transformations
- fit configured two-, three-, and four-peak candidates with deterministic multistart L-BFGS
- calculate BIC and AICc on the same deterministic grid of at most 600 points, using `4*k`
  parameters for `k` peaks
- retain the raw criterion minimum, but require convergence, full local Jacobian rank,
  acceptable column-normalized condition, no active bound, and at least 2% reconstructed area
  per component for structural selection
- integrate observed, reconstructed, and component curves by the trapezoidal rule on the full
  selected data interval
- estimate physical-parameter covariance from the numerical Jacobian and residual variance,
  and report two-sided Student-t intervals, correlation, rank, and condition diagnostics

These are local conditional intervals. Derivative residuals are serially correlated, and peak
count and preprocessing are selected rather than fixed by nature; the intervals therefore do
not represent unconditional coverage. BIC/AICc are comparison diagnostics under a common
sampling grid, not proof that each selected peak is a distinct reaction.

The evaluated joint model shares ordered-component skew across heating rates while keeping
height, center, and width experiment-specific. Its fit is compared with independent
three-peak fits on identical points. The real-data M5 audit rejects it as the default because
the predeclared convergence and 5% objective-increase gates are not met for every composition.
Input experiments are nevertheless canonically ordered, and synthetic permutation tests
establish numerical order invariance.

Legacy mapping: `fs_function.m`, `fs_mixture*.m`, `fs_fit.m`, `deconvolution*.m`, and `deconvolve*.m`.

## 7. Reaction-model layer

M6 estimates one overall-conversion kinetic triplet jointly across the four dynamic heating
programs for each dry-mass composition. It does not assign a common identity to any M5 peak.
The fit interval is the exact first-upward `alpha = 0.10` to `alpha = 0.90` crossing interval,
and each experiment contributes at most 250 deterministic points. Nonpositive measured rates
are counted and excluded because their logarithm is undefined; they are not clipped.

The fixed differential-model registry is:

| Name | `f(alpha)` |
|---|---|
| `f1`–`f4` | `(1-alpha)^q`, for `q = 1, 2, 3, 4` |
| `a2`–`a4` | `q*(1-alpha)*[-ln(1-alpha)]^((q-1)/q)`, for `q = 2, 3, 4` |
| `r2` | `2*(1-alpha)^(1/2)` |
| `r3` | `3*(1-alpha)^(2/3)` |
| `d1` | `1/(2*alpha)` |
| `d2` | `1/[-ln(1-alpha)]` |
| `d3` | `(3/2)*(1-alpha)^(2/3)/[1-(1-alpha)^(1/3)]` |
| `d4` | `(3/2)/[(1-alpha)^(-1/3)-1]` |
| `p2`–`p4` | `q*alpha^((q-1)/q)`, for `q = 2, 3, 4` |
| `random_scission` | `2*(sqrt(alpha)-alpha)` |

The simpler Šesták–Berggren form explored in the final scripts is named
`sestak_berggren_2`:

```text
f(alpha) = alpha^m * (1-alpha)^n
```

The generalized `sestak_berggren_3` form also includes

```text
f(alpha) = alpha^m * (1-alpha)^n * [-ln(1-alpha)]^p
```

Every reaction function requires `0 < alpha < 1`, is dimensionless, and must be positive and
finite. For fixed models, the regression is

```text
ln(dalpha/dt) - ln[f(alpha)] = ln(A/min^-1) - E/(R*T)
```

The empirical designs add the columns `ln(alpha)`, `ln(1-alpha)`, and optionally
`ln[-ln(1-alpha)]`. There is exactly one intercept, so `ln(A)` is identifiable rather than
being duplicated by an arbitrary amplitude coefficient. Least squares uses a QR-based solve,
not normal equations.

M6 retains all 19 candidates. BIC is the default comparison criterion. Full design rank, a
column-normalized condition threshold, `20 <= E <= 500 kJ/mol`, and `-5 <= m,n,p <= 5`
determine structural eligibility. The raw criterion minimum remains visible if a constraint
changes the selection. A criterion difference below two is flagged as ambiguous.

Each result reports rate and log-rate residuals, R², normalized rate RMSE, Durbin–Watson,
AICc/BIC, coefficient covariance and correlation, and local two-sided Student-t intervals.
Those intervals condition on the selected equation and an independent-error approximation;
they do not cover serial correlation, preprocessing, or model selection. Leave-one-heating-
rate-out fits therefore provide the primary M6 predictive gate, with accepted aggregate
log-rate RMSE no greater than `0.25`.

The compensation diagnostic fits every fixed model to a single heating-rate run, regresses
the resulting `ln(A)` values on `E`, and evaluates that line at an independently supplied
Starink reference energy. Its interval is a mean-response interval for that line. Because the
parameter pairs are estimated from the same data and kinetic compensation can have a
mathematical or error-induced origin, high R² is not treated as physical evidence and the
estimated `A` never replaces the direct multi-rate value.

The real-data audit selects empirical or third-order baselines depending on composition, but
all five selections fail the predeclared held-out-rate gate. A constant kinetic triplet is
therefore rejected as the default predictive model. M7 will compare these descriptive
baselines against a conversion-dependent isoconversional formulation.

Legacy mapping: `sb.m`, `sblin.m`, `sbrlin.m`, `comeffect*.m`, `compeffect*.m`, and `preexp*.m`.

## 8. Simulation contract

For a temperature program `T(t)`:

```text
dalpha/dt = A(alpha) * exp[-E(alpha)/(R*T(t))] * f(alpha)
```

The program may be a linear ramp, isothermal hold, or ramp-plus-hold. Interpolation of conversion-dependent kinetic parameters must define its support and boundary behavior. Solvers must preserve or explicitly terminate on the physical domain; taking `real(...)` of an invalid state is prohibited.

The 15 runs in `code/isotherm.mat` are ramp-and-hold validation experiments, not fitting data. Parameters estimated from the four dynamic heating-rate experiments must predict these programs without refitting. Required outputs are experimental/predicted conversion versus time, the measured and modeled temperature program, residuals separated into ramp and hold phases, and summary errors by composition and hold temperature.

Legacy mapping: `simulacionvyazov.m`, `simuisotherm.m`, `simsigmoid*.m`, and `movietga.m`.

## 9. Mixture-additivity and synergy contract

For waste-tire dry-mass fraction `x` and *Acrocomia aculeata* endocarp dry-mass fraction `1-x`, the legacy additive reference was

```text
alpha_additive(T) = (1-x)*alpha_endocarp(T) + x*alpha_tire(T)
```

with an analogous expression for conversion rate. The observed-minus-additive difference is a descriptive interaction metric. The new implementation must distinguish remaining mass, conversion, and rate metrics and propagate uncertainty in the component curves and composition.

Legacy mapping: `synergy.m`.

## 10. Primary references

- Vyazovkin, S. et al. (2011), “ICTAC Kinetics Committee recommendations for performing kinetic computations on thermal analysis data,” *Thermochimica Acta* 520, 1–19. <https://doi.org/10.1016/j.tca.2011.03.034>
- Friedman, H. L. (1964), “Kinetics of thermal degradation of char-forming plastics from thermogravimetry,” *Journal of Polymer Science Part C* 6, 183–195. <https://doi.org/10.1002/polc.5070060121>
- Ozawa, T. (1965), “A new method of analyzing thermogravimetric data,” *Bulletin of the Chemical Society of Japan* 38, 1881–1886. <https://doi.org/10.1246/bcsj.38.1881>
- Flynn, J. H. and Wall, L. A. (1966), “A quick, direct method for the determination of activation energy from thermogravimetric data,” *Journal of Polymer Science Part B* 4, 323–328. <https://doi.org/10.1002/pol.1966.110040504>
- Starink, M. J. (2003), “The determination of activation energy from linear heating rate experiments,” *Thermochimica Acta* 404, 163–176. <https://doi.org/10.1016/S0040-6031(03)00144-8>
- Vyazovkin, S. (2001), “Modification of the integral isoconversional method to account for variation in the activation energy,” *Journal of Computational Chemistry* 22, 178–183. <https://doi.org/10.1002/1096-987X(20010130)22:2%3C178::AID-JCC5%3E3.0.CO;2-%23>
- Vyazovkin, S. and Wight, C. A. (2000), “Estimating realistic confidence intervals for the activation energy determined from thermoanalytical measurements,” *Analytical Chemistry* 72, 3171–3175. <https://doi.org/10.1021/ac000210u>
- Fraser, R. D. B. and Suzuki, E. (1969), “Resolution of overlapping bands: Functions for simulating band shapes,” *Analytical Chemistry* 41, 37–39. <https://doi.org/10.1021/ac60270a007>
- Šesták, J. and Berggren, G. (1971), “Study of the kinetics of the mechanism of solid-state reactions at increasing temperatures,” *Thermochimica Acta* 3, 1–12. <https://doi.org/10.1016/0040-6031(71)85051-7>
- Roduit, B. (2000), “Computational aspects of kinetic analysis. Part E: The ICTAC Kinetics Project—numerical techniques and kinetics of solid state processes,” *Thermochimica Acta* 355, 171–180. <https://doi.org/10.1016/S0040-6031(00)00447-0>
- Vyazovkin, S. (2021), “A time to search: finding the meaning of the kinetic compensation effect,” *Molecules* 26, 3077. <https://doi.org/10.3390/molecules26113077>
- Barrie, P. J. (2012), “The mathematical origins of the kinetic compensation effect: 1. the effect of random experimental errors,” *Physical Chemistry Chemical Physics* 14, 318–326. <https://doi.org/10.1039/C1CP22666E>

## 11. Open method decisions

- The 150–700 °C coordinate is primary; 120 °C initial and 650/750/highest-common final
  profiles are required sensitivity results and are not folded into statistical intervals.
- FWO remains a compatibility method for thesis and literature comparison. Starink is the
  preferred linear integral approximation.
- Advanced-Vyazovkin uncertainty uses the pairwise-variance Fisher construction above. Its
  variance minimum is retained separately from the pairwise-objective energy estimate.
- Constant-triplet candidates remain M7 descriptive baselines because none passes the M6
  held-out-rate gate; evaluate conversion-dependent `E(alpha)` and `A(alpha)` for prediction.
- Decide whether DSC data are retained for diagnostics or excluded from version 1.0 analysis.
