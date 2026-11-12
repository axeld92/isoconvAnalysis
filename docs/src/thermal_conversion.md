# TG conversion and kinetics

Start with a [delimited thermal run](data_ingestion.md) that retains the measured signal
and the recorded temperature history. Select a program segment and the rows that support
the observation model. All row ranges below are inclusive and local to that segment.

## TG mass-loss conversion

For a remaining-mass signal `m`, conversion over a declared mass-loss event is
`α = (m_initial - m) / (m_initial - m_terminal)`. Use initial and terminal references
in the source unit. Reference windows estimate these values from the medians of the
recorded masses. The interval over which you evaluate conversion can differ from the
reference windows.

```julia
tg = tg_conversion(run, :tg;
    segment=1, reaction_rows=301:1500,
    initial_reference=251:300, terminal_reference=1501:1550)
```

Compare several plausible reference windows and reaction intervals before estimating
kinetics. A finite-temperature terminal mass may depend on the program and on the
instrument baseline. An apparent conversion outside `[0, 1]`, a local reversal, or a
missing crossing is informative. The conversion calculation preserves the recorded row
order and retains those observations. A signed change in mass divided by the original
sample mass gives a fractional mass loss on that sample basis. A reaction conversion
additionally requires a reference for the total mass lost in the event of interest.

### Common conversion and temperature support

Use `tg_conversion_crossings(tg, 0.5)` to inspect a target in one run. The diagnostic
follows acquisition order and reports upward crossings, downward crossings, isolated
touches, exact plateaus, and endpoint events. Linear interpolation within an observed
bracket uses the measured temperature and, where available, time. A plateau retains its
recorded interval and source lines.

To check a same-material group, prepare one `TGConversionResult` per run using a common
definition of the reaction and reference policy. Pass the results as `curves` with an
explicit conversion grid.

```julia
support = tg_conversion_support(curves, collect(0.10:0.05:0.90))
```

A target has common support when every run has one interior upward crossing and no
other event there. The group reports the minimum, mean, maximum, and span of the
interpolated crossing temperatures in kelvin. It also reports the reciprocal-temperature
span and the centered sum of squares of `1/T` for the contributing runs. A narrower span
or lower sum of squares can reduce regression leverage. The temperature statistics alone
leave the reaction model and each method's program assumptions to be checked separately.
Repeat the check with plausible reference windows and relevant program subsets.

### Method checks on TG data

At a unique upward crossing, `tg_friedman_rate` fits a local polynomial to the
recorded conversion against time. Declare a half-window in minutes and compare
several widths and degrees. A positive derivative allows `log(dα/dt)` to be
evaluated. Check how the differences between runs in that log rate vary across
the derivative settings. [ICTAC kinetic computation recommendations](https://doi.org/10.1016/j.tca.2011.03.034)

```julia
rate = tg_friedman_rate(tg, 0.50; half_window_min=0.4, degree=2)
```

For advanced Vyazovkin, use `tg_recorded_time_interval` to check both crossings and
the entire conversion branch between them. Strict support requires increasing
conversion throughout the selected interval. The diagnostic retains measured
interior time and temperature rows and identifies reversals and plateaus. Use
`tg_recorded_time_quadrature` on a supported interval to compare measured-node,
refined measured-time, and uniform-conversion quadrature at a declared energy probe.
Refinement tests the numerical integration of the piecewise-linear temperature
history. [Recorded-time method](https://doi.org/10.1002/%28SICI%291096-987X%28199702%2918%3A3%3C393%3A%3AAID-JCC9%3E3.0.CO%3B2-P) · [Incremental modification](https://doi.org/10.1002/1096-987X%2820010130%2922%3A2%3C178%3A%3AAID-JCC5%3E3.0.CO%3B2-%23)

```julia
interval = tg_recorded_time_interval(tg, 0.50; delta_alpha=0.02)
quadrature = tg_recorded_time_quadrature(interval, 120.0)
```

For KAS, FWO, and Starink, check the temperature history from a declared starting row
to the target crossing. `tg_ramp_diagnostic` compares the measured history with an
affine ramp over the same duration and endpoints. The probe output separates the
effect of the history, the temperature integral below the declared start, and the
shape error of each method's temperature-integral approximation. The energy probes
are inputs to these checks. [Comparison of linear-heating approximations](https://eprints.soton.ac.uk/18822/)

```julia
ramp = tg_ramp_diagnostic(tg, 0.50; start_row=first(tg.reaction_rows),
    probe_energies_kJ_per_mol=[60.0, 120.0, 240.0])
```

### Apparent activation energy

`tg_activation_energy` takes the checked TG curves and an explicit conversion grid.
Provide each run's local derivative half-window for Friedman or the numerical
configuration for advanced Vyazovkin. For the ramp comparisons, supply the declared
rates in the same order as `curves`.

```julia
targets = collect(0.10:0.10:0.90)
friedman = tg_activation_energy(:friedman, curves, targets;
    friedman_half_windows_min=half_windows_min)
settings = load_config("config/analysis_defaults.toml").isoconversional.advanced_vyazovkin
vyazovkin = tg_activation_energy(:advanced_vyazovkin, curves, targets;
    advanced=settings, time_subdivisions=8)
```

Each point retains its observed crossing support, run-specific statuses, and
regression or objective diagnostics. Inspect `status` before using an energy. A
blocked crossing, derivative, or interval leaves the energy absent. The numerical
search records flat, boundary, and multiple-minimum solutions for review. Compare
the primary settings with plausible alternatives on common support and report the
runs used. The public audit exercises these estimators on independent controls
before applying them to polypropylene.
