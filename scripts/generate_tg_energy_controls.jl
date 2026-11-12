module TGEnergyControls

using Printf
using SHA

const GAS_KJ_PER_MOL_K = 0.00831446261815324

# Integrate y = -log(1 - alpha) with RK4, separately from the kinetics package.
function integrated_conversion(times, derivative, subdivisions)
    hazard = 0.0
    conversion = zeros(length(times))
    for i in 2:length(times)
        start = times[i - 1]
        step = (times[i] - start) / subdivisions
        for j in 0:(subdivisions - 1)
            time = start + j * step
            k1 = derivative(time, hazard)
            k2 = derivative(time + step / 2, hazard + step * k1 / 2)
            k3 = derivative(time + step / 2, hazard + step * k2 / 2)
            k4 = derivative(time + step, hazard + step * k3)
            hazard += step * (k1 + 2k2 + 2k3 + k4) / 6
        end
        conversion[i] = -expm1(-hazard)
    end
    return conversion
end

function write_control(directory, protocol, control, rate, index)
    cold, hot = protocol["control_temperature_K"]
    duration = (hot - cold) / rate
    curvature = control["curved"] ? 0.2 * (-1)^(index - 1) : 0.0
    intercept = control["energy_intercept_kJ_per_mol"]
    slope = control["energy_slope_kJ_per_mol"]
    prefactor = protocol["control_prefactor_min_inv"]
    rows = protocol["control_rows"]
    initial_mass, terminal_mass = protocol["control_mass_mg"]

    function temperature(time)
        phase = time / duration
        return cold + (hot - cold) * (phase + curvature * sin(2pi * phase) / (2pi))
    end
    function derivative(time, hazard)
        energy = intercept - slope * expm1(-hazard)
        return prefactor * exp(-energy / (GAS_KJ_PER_MOL_K * temperature(time)))
    end

    times = [duration * (i / (rows - 1))^1.18 for i in 0:(rows - 1)]
    coarse = integrated_conversion(times, derivative, 4)
    fine = integrated_conversion(times, derivative, 8)
    difference = maximum(abs, fine .- coarse)
    material_id = control["id"]
    id = @sprintf("%s_%g", material_id, rate)
    difference <= protocol["control_refinement_tolerance_alpha"] &&
    last(fine) > maximum(protocol["alpha_grid"]) ||
        error("Control $id fails trajectory checks")

    filename = "$id.csv"
    path = joinpath(directory, filename)
    open(path, "w") do io
        println(io, "time_min,temperature_K,mass_mg")
        for (time, alpha) in zip(times, fine)
            mass = initial_mass - (initial_mass - terminal_mass) * alpha
            println(io, "$time,$(temperature(time)),$mass")
        end
    end
    return Dict{String,Any}(
        "id" => id,
        "material_id" => material_id,
        "filename" => filename,
        "sha256" => bytes2hex(open(sha256, path)),
        "rate_K_per_min" => rate,
        "curvature" => curvature,
        "refinement_max_abs_alpha" => difference,
        "last_alpha" => last(fine),
    )
end

function generate_tg_energy_controls(directory, protocol)
    runs = [
        write_control(directory, protocol, control, rate, index) for
        control in protocol["controls"] for
        (index, rate) in enumerate(protocol["control_rates_K_per_min"])
    ]
    println("Generated $(length(runs)) independently integrated TG controls.")
    return runs
end

end
