using OrdinaryDiffEqTsit5: Tsit5, solve
using SciMLBase: ODEProblem, successful_retcode

abstract type AbstractTemperatureProgram end

function _finite_value(value::Real, name::AbstractString)
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError("$name must be finite"))
    return converted
end

function _validated_program_time(program::AbstractTemperatureProgram, time_min::Real)
    time = _finite_value(time_min, "program time")
    duration = program_duration(program)
    time >= 0 && isapprox(time, 0.0; atol=1.0e-12, rtol=0.0) && return 0.0
    time >= 0 && isapprox(time, duration; atol=1.0e-12, rtol=1.0e-12) && return duration
    0 <= time <= duration ||
        throw(DomainError(time, "program time must lie in [0, $duration] min"))
    return time
end

"""
    SampledTemperatureProgram(time_min, temperature_K)

Continuous piecewise-linear temperature history, including cooling. Times must start at
zero and increase strictly. Temperatures are positive kelvin. The constructor copies both
input vectors.
"""
struct SampledTemperatureProgram <: AbstractTemperatureProgram
    time_min::Vector{Float64}
    temperature_K::Vector{Float64}

    function SampledTemperatureProgram(
        time_min::AbstractVector, temperature_K::AbstractVector
    )
        times = Float64.(time_min)
        values = Float64.(temperature_K)
        length(times) == length(values) > 0 ||
            throw(ArgumentError("temperature history must have equal nonempty channels"))
        all(isfinite, times) && first(times) == 0 && all(>(0), diff(times)) || throw(
            ArgumentError("history times must be finite, start at zero, and increase")
        )
        all(value -> isfinite(value) && value > 0, values) ||
            throw(ArgumentError("history temperatures must be finite and positive kelvin"))
        return new(times, values)
    end
end

"""
    program_duration(program)

Return the program duration in minutes.
"""
program_duration(program::SampledTemperatureProgram) = last(program.time_min)

"""
    temperature(program, time_min)

Evaluate the program at elapsed time `time_min` and return temperature in kelvin. Times
outside the program interval throw `DomainError`.
"""
function temperature(program::SampledTemperatureProgram, time_min::Real)
    time = _validated_program_time(program, time_min)
    index = searchsortedlast(program.time_min, time)
    index == length(program.time_min) && return last(program.temperature_K)
    fraction =
        (time - program.time_min[index]) /
        (program.time_min[index + 1] - program.time_min[index])
    return (1 - fraction) * program.temperature_K[index] +
           fraction * program.temperature_K[index + 1]
end

function _integrate_program_dose(rhs, sample_times; reltol, abstol)
    problem = ODEProblem(rhs, 0.0, (0.0, last(sample_times)))
    result = solve(
        problem,
        Tsit5();
        reltol,
        abstol,
        saveat=sample_times,
        save_everystep=false,
        dense=false,
    )
    successful_retcode(result) && length(result.u) == length(sample_times) ||
        error("recorded-time integration did not reach every requested observation")
    return result.u
end
