"""
Fixed-energy integration diagnostics for one recorded TG program prefix. The
measured and affine integrals use time in minutes and the same start and end
temperatures. `omitted_start_fraction` is the fraction of the affine
zero-temperature integral below the selected start temperature.
`approximation_log_residuals` contains the KAS, FWO and Starink log-shape
residuals relative to the numerical temperature integral. Temperature-independent
constants in those approximations are omitted. Differences across runs at the
same probe energy cancel those constants.
"""
struct TGRampIntegralProbe
    energy_kJ_per_mol::Float64
    crossing_E_over_RT::Float64
    log_measured_integral_min::Float64
    log_affine_integral_min::Float64
    measured_vs_affine_log_ratio::Float64
    max_quadrature_log_difference::Float64
    omitted_start_fraction::Float64
    approximation_log_residuals::NamedTuple{(:kas, :fwo, :starink),NTuple{3,Float64}}
end

"""
Diagnostics for the history from a declared TG row to an observed target crossing.
`status` is `:evaluated` for a timed prefix with positive net heating, or identifies
why it could not be evaluated. The sample-weighted OLS slope and the endpoint
secant describe the observed history. Residuals use the endpoint affine control.
The half-duration secants split the history by elapsed time. No status imposes a
numerical tolerance on ramp departure.
"""
struct TGRampDiagnostic
    alpha::Float64
    status::Symbol
    crossing::TGCrossingDiagnostic
    start_row::Int
    start_source_line::Int
    start_time_min::Union{Nothing,Float64}
    start_temperature_K::Float64
    duration_min::Union{Nothing,Float64}
    secant_K_per_min::Union{Nothing,Float64}
    OLS_K_per_min::Union{Nothing,Float64}
    affine_max_abs_residual_K::Union{Nothing,Float64}
    affine_RMSE_K::Union{Nothing,Float64}
    half_secants_K_per_min::Union{Nothing,Tuple{Float64,Float64}}
    step_slope_range_K_per_min::Union{Nothing,Tuple{Float64,Float64}}
    falling_temperature_steps::Union{Nothing,Int}
    flat_temperature_steps::Union{Nothing,Int}
    cumulative_cooling_K::Union{Nothing,Float64}
    raw_row_count::Int
    probes::Vector{TGRampIntegralProbe}
end

# Positive abscissae and corresponding weights for symmetric Gauss–Legendre rules.
const _RAMP_GAUSS_4 = (
    (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374539)
)
const _RAMP_GAUSS_8 = (
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
)

function _ramp_gauss(f, lo, hi, rule)
    half_width = (hi - lo) / 2
    center = lo + half_width
    return half_width * sum(rule) do (node, weight)
        displacement = half_width * node
        return weight * (f(center - displacement) + f(center + displacement))
    end
end

function _ramp_log_integral(t, T, scaled_energy, shift, rule)
    kernel(temperature) = exp(-scaled_energy / temperature - shift)
    value = sum(2:length(t)) do i
        lo, hi = minmax(T[i - 1], T[i])
        if lo == hi
            return (t[i] - t[i - 1]) * kernel(lo)
        end
        return (t[i] - t[i - 1]) * (_ramp_gauss(kernel, lo, hi, rule) / (hi - lo))
    end
    isfinite(value) && value > 0 ||
        throw(ArgumentError("ramp probe integral is outside the representable range"))
    return shift + log(value)
end

# p(x) = integral from x to infinity of exp(-u)/u^2 du.
function _ramp_log_temperature_integral(x)
    isfinite(x) && x > 0 ||
        throw(ArgumentError("ramp probe E/(RT) must be positive and finite"))
    if x >= 1
        kernel(s) = exp(-s) * (x / (x + s))^2
        # After scaling by exp(-x)/x^2, the omitted tail is bounded by exp(-48).
        scaled = sum(_ramp_gauss(kernel, i, i + 1, _RAMP_GAUSS_8) for i in 0:47)
        return -x - 2log(x) + log(scaled)
    end
    # Substitute u = x/T_fraction. Geometric subintervals resolve the boundary
    # layer near zero when x is small.
    kernel_fraction(y) = exp(-x / y + x)
    integrate_fraction(a, b) = begin
        edges = range(a, b; length=9)
        sum(
            _ramp_gauss(kernel_fraction, edges[i], edges[i + 1], _RAMP_GAUSS_8) for
            i in 1:8
        )
    end
    lo = max(x, eps(Float64))
    scaled = integrate_fraction(0.0, lo)
    while lo < 1
        hi = min(2lo, 1.0)
        scaled += integrate_fraction(lo, hi)
        lo = hi
    end
    return -x - log(x) + log(scaled)
end

function _ramp_integral_probe(energy, t, T, affine)
    θ = energy / GAS_CONSTANT_KJ_PER_MOL_K
    shift = -θ / maximum(T)
    measured = _ramp_log_integral(t, T, θ, shift, _RAMP_GAUSS_8)
    control = _ramp_log_integral(t, affine, θ, shift, _RAMP_GAUSS_8)
    measured_coarse = _ramp_log_integral(t, T, θ, shift, _RAMP_GAUSS_4)
    control_coarse = _ramp_log_integral(t, affine, θ, shift, _RAMP_GAUSS_4)
    x = θ / last(T)
    log_exact = _ramp_log_temperature_integral(x)
    log_start = _ramp_log_temperature_integral(θ / first(T))
    return TGRampIntegralProbe(
        energy,
        x,
        measured,
        control,
        measured - control,
        max(abs(measured - measured_coarse), abs(control - control_coarse)),
        exp(log_start - log_exact),
        (
            kas=(-x - 2log(x) - log_exact),
            fwo=-1.052x - log_exact,
            starink=-1.0008x - 1.92log(x) - log_exact,
        ),
    )
end

"""
    tg_ramp_diagnostic(curve, alpha; start_row, probe_energies_kJ_per_mol=Float64[])

Check the measured temperature program from an explicit segment-local start row
through the unique upward crossing at `alpha`. Use the recorded clock and
piecewise-linear temperature history, including the interpolated crossing
endpoint. Compare the measured prefix with an affine ramp through the same
endpoints and elapsed time. Optional finite positive probe energies quantify
Arrhenius integral sensitivity and separate the formula shapes from the program
departure. The probes are specified inputs, with no kinetic parameter fitted.
"""
function tg_ramp_diagnostic(
    curve::TGConversionResult,
    alpha::Real;
    start_row::Integer,
    probe_energies_kJ_per_mol::AbstractVector{<:Real}=Float64[],
)
    diagnostic = tg_conversion_crossings(curve, alpha)
    start_row in curve.reaction_rows ||
        throw(ArgumentError("ramp start row must be within the selected TG interval"))
    energies = Float64.(probe_energies_kJ_per_mol)
    all(e -> isfinite(e) && e > 0, energies) && allunique(energies) ||
        throw(ArgumentError("ramp probe energies must be finite, positive and distinct"))
    clock = curve.time_min
    offset = first(curve.reaction_rows) - 1
    first_index = start_row - offset
    start_time = clock === nothing ? nothing : clock[first_index]
    start_T = curve.temperature_K[first_index]
    prefix = (
        diagnostic.alpha,
        diagnostic,
        Int(start_row),
        curve.source_lines[first_index],
        start_time,
        start_T,
    )
    function unavailable(status)
        target, crossing, start, source, time, temperature = prefix
        return TGRampDiagnostic(
            target,
            status,
            crossing,
            start,
            source,
            time,
            temperature,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            0,
            TGRampIntegralProbe[],
        )
    end
    diagnostic.status === :unique_upward || return unavailable(:unsupported_crossing)
    clock === nothing && return unavailable(:missing_time)
    event = only(diagnostic.crossings)
    crossing_time = event.time_min::Float64
    crossing_time > start_time || return unavailable(:start_at_or_after_crossing)
    last_index = searchsortedlast(clock, crossing_time)
    indices = first_index:last_index
    t = clock[indices] .- start_time
    T = curve.temperature_K[indices]
    duration = crossing_time - start_time
    if clock[last_index] < crossing_time
        push!(t, duration)
        push!(T, event.temperature_K::Float64)
    else
        T[end] = event.temperature_K::Float64
    end
    dT = diff(T)
    secant = (last(T) - first(T)) / duration
    centered_t = t .- mean(t)
    ols = sum(centered_t .* (T .- mean(T))) / sum(abs2, centered_t)
    affine = first(T) .+ secant .* t
    affine[end] = last(T)
    departure = T .- affine
    half_time = duration / 2
    mid_index = searchsortedlast(t, half_time)
    mid_fraction = (half_time - t[mid_index]) / (t[mid_index + 1] - t[mid_index])
    mid_T = T[mid_index] + mid_fraction * (T[mid_index + 1] - T[mid_index])
    half_secants = ((mid_T - first(T)) / half_time, (last(T) - mid_T) / half_time)
    probes = if secant > 0
        [_ramp_integral_probe(energy, t, T, affine) for energy in energies]
    else
        TGRampIntegralProbe[]
    end
    target, crossing, start, source, time, temperature = prefix
    return TGRampDiagnostic(
        target,
        secant > 0 ? :evaluated : :nonpositive_net_ramp,
        crossing,
        start,
        source,
        time,
        temperature,
        duration,
        secant,
        ols,
        maximum(abs, departure),
        sqrt(mean(abs2, departure)),
        half_secants,
        extrema(dT ./ diff(t)),
        count(<(0), dT),
        count(iszero, dT),
        sum(x -> max(0.0, -x), dT),
        length(indices),
        probes,
    )
end
