"""
Local conversion-rate evidence at one observed TG crossing. `status` is
`:positive_rate`, `:nonpositive_rate`, `:crossing_unavailable`, `:time_unavailable`,
`:insufficient_window`, or `:fit_unavailable`. `crossing_status` gives the raw
crossing classification. A fitted rate is signed and expressed per minute.
`window_source_lines` bounds the physical source lines used for the fit, and
`residual_rmse` describes the local conversion fit without a noise model.
"""
struct TGFriedmanRateDiagnostic
    alpha::Float64
    crossing_status::Symbol
    status::Symbol
    crossing_temperature_K::Union{Nothing,Float64}
    crossing_time_min::Union{Nothing,Float64}
    crossing_source_lines::Union{Nothing,Tuple{Int,Int}}
    half_window_min::Float64
    degree::Int
    window_source_lines::Union{Nothing,Tuple{Int,Int}}
    window_time_min::Union{Nothing,Tuple{Float64,Float64}}
    window_conversion_range::Union{Nothing,Tuple{Float64,Float64}}
    points_before::Int
    points_after::Int
    window_points::Int
    conversion_rate_min_inv::Union{Nothing,Float64}
    residual_rmse::Union{Nothing,Float64}
end

"""
    tg_friedman_rate(curve, alpha; half_window_min, degree=2, minimum_points_per_side=2)

Estimate local `dα/dt` at the time of a unique upward crossing in `curve`. Use
recorded conversion values from a fixed time window centered at the interpolated
crossing time. The local polynomial is fitted in scaled time with Gaussian
weights `exp(-0.5 * ((t - t_cross) / half_window_min)^2)`. Degrees 1 through 3
are supported. The fit requires at least `degree + 2` total samples and the
requested minimum number of observations strictly on each side of the crossing.
The requested window is retained without expansion.

The reported rate retains its sign. A positive result makes the local logarithm
defined for these settings. Compare settings and measurement support separately
before using it in a Friedman regression.
"""
function tg_friedman_rate(
    curve::TGConversionResult,
    alpha::Real;
    half_window_min::Real,
    degree::Integer=2,
    minimum_points_per_side::Integer=2,
)
    width = Float64(half_window_min)
    isfinite(width) && width > 0 ||
        throw(ArgumentError("half_window_min must be finite and positive"))
    1 <= degree <= 3 ||
        throw(ArgumentError("Friedman local polynomial degree must be in 1:3"))
    minimum_points_per_side >= 1 ||
        throw(ArgumentError("minimum_points_per_side must be positive"))

    diagnostic = tg_conversion_crossings(curve, alpha)
    crossing = diagnostic.status === :unique_upward ? only(diagnostic.crossings) : nothing
    T = crossing === nothing ? nothing : crossing.temperature_K
    t = crossing === nothing ? nothing : crossing.time_min
    lines = crossing === nothing ? nothing : crossing.source_lines
    result(
        status;
        window=nothing,
        time_range=nothing,
        conversion_range=nothing,
        before=0,
        after=0,
        n=0,
        rate=nothing,
        rmse=nothing,
    ) = TGFriedmanRateDiagnostic(
        diagnostic.alpha,
        diagnostic.status,
        status,
        T,
        t,
        lines,
        width,
        Int(degree),
        window,
        time_range,
        conversion_range,
        before,
        after,
        n,
        rate,
        rmse,
    )
    crossing === nothing && return result(:crossing_unavailable)
    t === nothing && return result(:time_unavailable)

    time = curve.time_min::Vector{Float64}
    rows = searchsortedfirst(time, t - width):searchsortedlast(time, t + width)
    isempty(rows) && return result(:insufficient_window)
    window = (curve.source_lines[first(rows)], curve.source_lines[last(rows)])
    scaled_time = (view(time, rows) .- t) ./ width
    before = count(<(0), scaled_time)
    after = count(>(0), scaled_time)
    n = length(rows)
    time_range = (time[first(rows)], time[last(rows)])
    conversion_range = extrema(view(curve.alpha, rows))
    provenance = (; window, time_range, conversion_range, before, after, n)
    if min(before, after) < minimum_points_per_side || n < degree + 2
        return result(:insufficient_window; provenance...)
    end

    response = view(curve.alpha, rows) .- diagnostic.alpha
    design = reduce(hcat, (scaled_time .^ power for power in 0:degree))
    sqrt_weight = @. exp(-0.25 * scaled_time^2)
    weighted_design = design .* sqrt_weight
    rank(weighted_design) == degree + 1 || return result(:fit_unavailable; provenance...)
    coefficients = weighted_design \ (response .* sqrt_weight)
    rate = coefficients[2] / width
    rmse = sqrt(mean(abs2, design * coefficients .- response))
    isfinite(rate) && isfinite(rmse) || return result(:fit_unavailable; provenance...)
    status = rate > 0 ? :positive_rate : :nonpositive_rate
    return result(status; provenance..., rate, rmse)
end
