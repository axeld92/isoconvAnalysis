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
