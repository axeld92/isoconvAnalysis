# Fraser–Suzuki deconvolution

M5 decomposes each processed dynamic `dalpha/dT` curve into an ordered sum of two, three,
or four Fraser–Suzuki peaks. The decomposition is descriptive: an ordered peak is not a
chemical species or elementary reaction unless independent evidence establishes that
interpretation.

## Peak convention

`FraserSuzukiPeak` stores height in `K^-1`, dimensionless skew, center temperature in kelvin,
and width in kelvin. Width is the Gaussian full width at half maximum when skew is zero. The
implementation evaluates the Gaussian limit near zero skew and returns the continuous value
zero outside the logarithm's real domain.

```julia
peak = FraserSuzukiPeak(0.004, -0.2, 620.0, 80.0)
rate = fraser_suzuki(collect(500.0:5.0:800.0), peak)
```

## Independent fits and peak-count selection

The accepted default fits the first upward conversion interval from `alpha = 0.05` to
`alpha = 0.95`. Heights are nonnegative, widths and skews are bounded, and centers are
ordered with a minimum separation. These constraints are encoded with smooth parameter
transformations; deterministic multistart L-BFGS optimization is used for a fixed input and
configuration.

```julia
config = load_config("config/analysis_defaults.toml")

fit = fit_deconvolution(processed_experiment, 3, config)
comparison = compare_peak_counts(processed_experiment, config)
selected = comparison.results[
    findfirst(result -> result.peak_count == comparison.selected_peak_count,
              comparison.results)
]
```

`compare_peak_counts` retains all candidate fits. It reports both the raw BIC or AICc minimum
and the accepted structural selection. A candidate is structurally ineligible if it did not
converge, its local Jacobian is rank or condition deficient, a parameter or center separation
is active at a configured bound, or a component contributes less than the configured minimum
area fraction. A difference below two criterion units is reported as ambiguous.
If every candidate is ineligible, the raw criterion minimum is retained with an explicit
warning so the experiment remains visible but is not misrepresented as structurally accepted.

## Result evidence

`DeconvolutionResult` contains the selected temperature interval, observed and reconstructed
rates, every component curve, integrated component areas and fractions, the physical-parameter
covariance and correlation matrices, local Student-t intervals, optimizer outcomes, residual
metrics, Durbin–Watson statistic, boundary flags, and a provenance fingerprint. Negative
measured derivative values remain in the residual calculation; the fitter never clips them.

The parameter intervals are local linearized intervals conditional on the selected peak
count and preprocessing. They do not account for residual autocorrelation, model-selection
uncertainty, or alternative preprocessing profiles. High correlations and active bounds are
therefore substantive warnings, not cosmetic optimizer messages.

## Joint multi-rate diagnostic

```julia
joint = fit_joint_deconvolution(composition_group, 3, config)
```

The joint formulation canonically orders experiments and shares skew by ordered component,
while height, center, and width remain rate-specific. M5 keeps this API for sensitivity work,
but the real-data gate did not support making it the default: two of five composition fits
did not converge, and two composition groups exceeded the predeclared 5% objective-increase
threshold relative to independent fits.

See the complete real-data audit in `docs/deconvolution_audit.md` and its machine-readable
evidence in `docs/audits/m5_deconvolution_results.toml`. The audit covers all 20 dynamic
calibration runs and leaves all 15 ramp-and-hold programs held out.

The exported types and methods are included in the manual's generated API index.
