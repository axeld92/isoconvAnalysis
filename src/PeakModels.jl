"""
One Fraser–Suzuki peak with height, dimensionless skew, center temperature in kelvin, and
width in kelvin. The width is the full width at half maximum in the Gaussian limit.
"""
struct FraserSuzukiPeak
    height_K_inv::Float64
    skew::Float64
    center_K::Float64
    width_K::Float64

    function FraserSuzukiPeak(height::Real, skew::Real, center::Real, width::Real)
        values = Float64.((height, skew, center, width))
        all(isfinite, values) || throw(ArgumentError("peak parameters must be finite"))
        values[1] >= 0 || throw(ArgumentError("peak height must be nonnegative"))
        values[4] > 0 || throw(ArgumentError("peak width must be positive"))
        return new(values...)
    end
end

function _fraser_suzuki_value(temperature, height, skew, center, width)
    if abs(skew) <= 1.0e-7
        coordinate = (temperature - center) / width
        return height * exp(-4 * log(2) * coordinate^2)
    end
    argument = 1 + 2 * skew * (temperature - center) / width
    argument > 0 || return zero(temperature + height + skew + center + width)
    exponent = -log(2) * (log(argument) / skew)^2
    return height * exp(exponent)
end

"""
    fraser_suzuki(temperature_K, peak)

Evaluate a Fraser–Suzuki peak. Values outside the logarithm's real domain are exactly zero,
which is the continuous limiting value at the domain boundary. The numerically stable
Gaussian limit is used near zero skew.
"""
function fraser_suzuki(temperature_K::Real, peak::FraserSuzukiPeak)
    return Float64(
        _fraser_suzuki_value(
            Float64(temperature_K),
            peak.height_K_inv,
            peak.skew,
            peak.center_K,
            peak.width_K,
        ),
    )
end

function fraser_suzuki(temperature_K::AbstractVector{<:Real}, peak::FraserSuzukiPeak)
    return [fraser_suzuki(temperature, peak) for temperature in temperature_K]
end

"""
    fraser_suzuki_mixture(temperature_K, peaks)

Evaluate the additive curve for an ordered collection of Fraser–Suzuki peaks.
"""
function fraser_suzuki_mixture(
    temperature_K::AbstractVector{<:Real}, peaks::AbstractVector{FraserSuzukiPeak}
)
    isempty(peaks) && throw(ArgumentError("a mixture requires at least one peak"))
    result = zeros(Float64, length(temperature_K))
    for peak in peaks
        result .+= fraser_suzuki(temperature_K, peak)
    end
    return result
end
