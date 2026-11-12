"""
One apparent TG activation-energy estimate with common crossing support.
`run_statuses` follows the result's run order. `method_inputs` contains the local
rates for Friedman, declared ramp rates for a ramp approximation, or log time
integrals at the optimum for advanced Vyazovkin. `status` exposes blocked
observations and numerical warnings. A blocked point has no energy. Numerical
warnings can accompany a reported estimate. `regression` and `objective` contain
descriptive fit diagnostics for the applicable method.
"""
struct TGActivationEnergyPoint
    support::TGConversionSupportPoint
    status::Symbol
    energy_kJ_per_mol::Union{Nothing,Float64}
    run_statuses::Vector{Symbol}
    method_inputs::Vector{Float64}
    regression::Union{Nothing,LinearRegressionDiagnostics}
    objective::Union{Nothing,Float64}
    objective_relative_range::Union{Nothing,Float64}
    objective_local_minima::Union{Nothing,Int}
end

function _tg_unavailable_energy(point, status, run_statuses)
    return TGActivationEnergyPoint(
        point, status, nothing, run_statuses, Float64[], nothing, nothing, nothing, nothing
    )
end

function _tg_positive_inputs(values, count, label)
    values === nothing && throw(ArgumentError("declare $label in curve order"))
    inputs = Float64.(collect(values))
    length(inputs) == count && all(x -> isfinite(x) && x > 0, inputs) ||
        throw(ArgumentError("$label must contain one positive finite value per curve"))
    return inputs
end

function _tg_advanced_settings(advanced, time_subdivisions)
    advanced isa AdvancedVyazovkinConfig ||
        throw(ArgumentError("declare AdvancedVyazovkinConfig for recorded-time estimation"))
    lower, upper = advanced.energy_bounds_kJ_per_mol
    0 < advanced.delta_alpha < 1 &&
    isfinite(lower) &&
    isfinite(upper) &&
    lower < upper &&
    advanced.optimizer_grid_points >= 3 &&
    isfinite(advanced.optimizer_tolerance_kJ_per_mol) &&
    advanced.optimizer_tolerance_kJ_per_mol > 0 &&
    advanced.optimizer_max_iterations > 0 &&
    isfinite(advanced.boundary_tolerance_kJ_per_mol) &&
    0 <= advanced.boundary_tolerance_kJ_per_mol < (upper - lower) / 2 &&
    isfinite(advanced.flat_relative_tolerance) &&
    advanced.flat_relative_tolerance >= 0 ||
        throw(ArgumentError("invalid recorded-time energy search settings"))
    time_subdivisions >= 1 || throw(ArgumentError("time_subdivisions must be positive"))
    return advanced
end

function _tg_advanced_energy(point, curves, advanced, subdivisions)
    intervals = [
        tg_recorded_time_interval(curve, point.alpha; delta_alpha=advanced.delta_alpha) for
        curve in curves
    ]
    statuses = getfield.(intervals, :status)
    all(==(:supported), statuses) ||
        return _tg_unavailable_energy(point, :unsupported_interval, statuses)

    log_integrals(energy) = [
        _tg_log_time_subdivided_integral(
            energy, interval.time_min, interval.temperature_K, subdivisions
        ) for interval in intervals
    ]
    objective(energy) = pairwise_vyazovkin_objective(log_integrals(energy))
    minimum = _scan_and_minimize(objective, advanced)
    lower, upper = advanced.energy_bounds_kJ_per_mol
    energy = minimum.energy
    status = if minimum.relative_range <= advanced.flat_relative_tolerance
        :flat_objective
    elseif min(energy - lower, upper - energy) <= advanced.boundary_tolerance_kJ_per_mol
        :boundary_solution
    elseif minimum.local_minima > 1
        :multiple_minima
    else
        :estimated
    end
    return TGActivationEnergyPoint(
        point,
        status,
        energy,
        statuses,
        log_integrals(energy),
        nothing,
        minimum.value,
        minimum.relative_range,
        minimum.local_minima,
    )
end

function _tg_linear_energy(method, point, curves, inputs, degree, minimum_points)
    temperatures = Float64[only(d.crossings).temperature_K for d in point.crossings]
    if method === :friedman
        rates = [
            tg_friedman_rate(
                curve,
                point.alpha;
                half_window_min=width,
                degree,
                minimum_points_per_side=minimum_points,
            ) for (curve, width) in zip(curves, inputs)
        ]
        statuses = getfield.(rates, :status)
        all(==(:positive_rate), statuses) ||
            return _tg_unavailable_energy(point, :unsupported_rate, statuses)
        values = Float64[r.conversion_rate_min_inv for r in rates]
        response, coefficient = log.(values), 1.0
    else
        statuses = fill(:declared_ramp_rate, length(curves))
        values = inputs
        response, coefficient = _ramp_method_response(method, temperatures, values)
    end
    regression = _linear_regression(inv.(temperatures), response)
    energy = -GAS_CONSTANT_KJ_PER_MOL_K * regression.slope / coefficient
    status = energy > 0 ? :estimated : :nonpositive_energy
    return TGActivationEnergyPoint(
        point, status, energy, statuses, copy(values), regression, nothing, nothing, nothing
    )
end

"""
    tg_activation_energy(method, curves, targets; material_id=nothing,
                         friedman_half_windows_min=nothing, friedman_degree=2,
                         minimum_points_per_side=2, advanced=nothing,
                         time_subdivisions=8, ramp_rates_K_per_min=nothing)

Estimate apparent `E(alpha)` in kJ/mol at interior targets for three or more
same-material TG runs. Returns the method, material, ordered run IDs and a vector
of [`TGActivationEnergyPoint`](@ref). A target requires a unique upward crossing
in every run. Linear estimators also require nonzero reciprocal-temperature spread.
All input vectors follow the order of `curves`.

Use `:friedman` with declared derivative half-windows. Use `:advanced_vyazovkin`
with `advanced` search and conversion-increment settings. The advanced objective
integrates the recorded, piecewise-linear time paths with `time_subdivisions`
per measured step. Use `:kas`, `:fwo`, or `:starink` with explicitly declared
constant-ramp rates after checking the approximation on the recorded program.
The fit diagnostics describe the fitted point and do not include an uncertainty
interval or establish a reaction mechanism.
"""
function tg_activation_energy(
    method::Symbol,
    curves::AbstractVector{TGConversionResult},
    targets::AbstractVector{<:Real};
    material_id=nothing,
    friedman_half_windows_min=nothing,
    friedman_degree::Integer=2,
    minimum_points_per_side::Integer=2,
    advanced=nothing,
    time_subdivisions::Integer=8,
    ramp_rates_K_per_min=nothing,
)
    method in _SUPPORTED_ISOCONVERSIONAL_METHODS ||
        throw(ArgumentError("unsupported TG method $method"))
    length(curves) >= 3 ||
        throw(ArgumentError("TG estimation requires at least three runs"))
    all(a -> 0 < a < 1, targets) ||
        throw(ArgumentError("energy targets must lie in (0, 1)"))
    support = tg_conversion_support(curves, targets; material_id)
    inputs = if method === :friedman
        1 <= friedman_degree <= 3 && minimum_points_per_side >= 1 ||
            throw(ArgumentError("invalid local derivative settings"))
        _tg_positive_inputs(
            friedman_half_windows_min, length(curves), "derivative half-windows"
        )
    elseif method === :advanced_vyazovkin
        _tg_advanced_settings(advanced, time_subdivisions)
    else
        _tg_positive_inputs(ramp_rates_K_per_min, length(curves), "constant-ramp rates")
    end
    points = map(support.points) do point
        statuses = getfield.(point.crossings, :status)
        if !point.common_upward
            _tg_unavailable_energy(point, :unsupported_crossing, statuses)
        elseif method === :advanced_vyazovkin
            _tg_advanced_energy(point, curves, inputs, time_subdivisions)
        elseif point.reciprocal_temperature_sxx_K_inv2 == 0
            _tg_unavailable_energy(point, :zero_temperature_spread, statuses)
        else
            _tg_linear_energy(
                method, point, curves, inputs, friedman_degree, minimum_points_per_side
            )
        end
    end
    return (; method, material_id=support.material_id, run_ids=support.run_ids, points)
end
