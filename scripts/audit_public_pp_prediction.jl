using IsoconversionalAnalysis
using IsoconversionalAnalysis:
    Optim, SampledTemperatureProgram, _integrate_program_dose, temperature
using SHA
using Statistics
using TOML

module PublicPPPredictionInputs
include("audit_public_tg_conversion.jl")
end
using .PublicPPPredictionInputs:
    TG_ROOT, TG_IMPORT_MANIFEST, TG_PROTOCOL, import_tg, tg_case_curves

const PP_PROTOCOL = "config/public_pp_prediction.toml"
const PP_REPORT = "docs/audits/public_pp_prediction.toml"
const PP_MODELS = (:constant_energy, :conversion_dependent_energy)
const PP_R = 0.00831446261815324

struct PPProgram
    history::SampledTemperatureProgram
    first_mass_mg::Float64
end

struct PPObservation
    id::String
    sha256::String
    program::PPProgram
    mass_mg::Vector{Float64}
    reaction_rows::UnitRange{Int}
    conversion::Union{Nothing,TGConversionResult}
end

function weights_for_time(times)
    steps = diff(times)
    length(steps) >= 1 && all(>(0), steps) || error("a weighted window needs ordered time")
    weights = [first(steps); steps[1:(end - 1)] .+ steps[2:end]; last(steps)]
    return weights ./ (2 * (last(times) - first(times)))
end

function fit_rows(time, rows, count)
    samples = view(time, rows)
    grid = range(first(samples), last(samples); length=min(length(rows), count))
    return unique!([first(rows); searchsortedfirst.(Ref(time), grid); last(rows)])
end

function calibration_record(observation, count)
    t = observation.program.history.time_min
    full = fit_rows(t, 1:length(t), count)
    reaction = fit_rows(t, observation.reaction_rows, count)
    rows = sort!(union(full, reaction))
    w = zeros(length(rows))
    w[searchsortedfirst.(Ref(rows), full)] .+= weights_for_time(t[full]) ./ 2
    w[searchsortedfirst.(Ref(rows), reaction)] .+= weights_for_time(t[reaction]) ./ 2
    return (;
        program=observation.program,
        t=t[rows],
        y=observation.mass_mg[rows] ./ observation.program.first_mass_mg,
        w,
    )
end

# Integrated dose h has dh/dt = k(α,T). The closed-form inverse of the integrated
# nth-order reaction function supplies remaining reactive mass without an ODE singularity.
function reactive_fraction(dose, n)
    h = max(0.0, dose)
    q = n - 1
    q == 0 && return exp(-h)
    q * h <= -1 && return 0.0
    return exp(-log1p(q * h) / q)
end

function energy_offset(shape, alpha)
    shape === nothing && return 0.0
    grid, offsets = shape
    alpha <= first(grid) && return first(offsets)
    alpha >= last(grid) && return last(offsets)
    i = searchsortedlast(grid, alpha)
    w = (alpha - grid[i]) / (grid[i + 1] - grid[i])
    return (1 - w) * offsets[i] + w * offsets[i + 1]
end

function reactive_history(
    program::PPProgram, sample_times, parameters, shape, protocol; refined=false
)
    E, log_rate, n = parameters
    reference = inv(protocol["reference_temperature_K"])
    history = program.history
    function rhs(h, _, time)
        offset = energy_offset(shape, 1 - reactive_fraction(h, n))
        return exp(
            log_rate - (E + offset) * (inv(temperature(history, time)) - reference) / PP_R
        )
    end
    reltol = protocol[refined ? "check_relative_tolerance" : "relative_tolerance"]
    abstol = protocol[refined ? "check_absolute_tolerance" : "absolute_tolerance"]
    dose = _integrate_program_dose(rhs, sample_times; reltol, abstol)
    return [reactive_fraction(h, n) for h in dose]
end

function residue_fit(pieces, active, upper)
    numerator, denominator = 0.0, 0.0
    for (piece, survival) in zip(pieces, active)
        alpha = 1 .- survival
        numerator += sum(piece.w .* alpha .* (piece.y .- survival))
        denominator += sum(piece.w .* abs2.(alpha))
    end
    return denominator == 0 ? 0.0 : clamp(numerator / denominator, 0.0, upper)
end

function pp_training_objective(pieces, parameters, shape, protocol; refined=false)
    reactive = [
        reactive_history(piece.program, piece.t, parameters, shape, protocol; refined) for
        piece in pieces
    ]
    residue = residue_fit(pieces, reactive, protocol["max_retained_fraction"])
    loss = mean(zip(pieces, reactive)) do (piece, survival)
        predicted = residue .+ (1 - residue) .* survival
        return sum(piece.w .* abs2.(predicted .- piece.y))
    end
    return loss, residue
end

function fit_pp(training, shape, protocol)
    pieces = [
        calibration_record(run, protocol["fit_points_per_window"]) for run in training
    ]
    bounds = [
        protocol["energy_bounds_kJ_per_mol"],
        protocol["log_reference_rate_bounds_min_inv"],
        protocol["reaction_order_bounds"],
    ]
    lower, upper = first.(bounds), last.(bounds)
    width = upper .- lower
    decode(x) = lower .+ x .* width
    function minimize(start; refined=false)
        function objective(x)
            all(y -> 0 <= y <= 1, x) || return 1e6
            return first(pp_training_objective(pieces, decode(x), shape, protocol; refined))
        end
        initial = (start .- lower) ./ width
        fit = optimize(
            objective,
            initial,
            NelderMead(;
                initial_simplex=Optim.AffineSimplexer(; a=refined ? 0.005 : 0.06, b=0.0)
            ),
            Optim.Options(;
                iterations=protocol[if refined
                    "refinement_iterations"
                else
                    "optimizer_iterations"
                end],
                g_abstol=protocol[if refined
                    "refinement_objective_tolerance"
                else
                    "optimizer_objective_tolerance"
                end],
            ),
        )
        return (;
            parameters=decode(Optim.minimizer(fit)),
            value=Optim.minimum(fit),
            converged=Optim.converged(fit),
            iterations=Optim.iterations(fit),
        )
    end
    results = minimize.(protocol["starts"])
    usable = filter(fit -> fit.converged && fit.value < 1e6, results)
    isempty(usable) && error("no converged mass fit for $(getfield.(training, :id))")
    best = usable[argmin(getfield.(usable, :value))]
    refinement = minimize(best.parameters; refined=true)
    baseline = pp_training_objective(pieces, best.parameters, shape, protocol; refined=true)
    parameters =
        refinement.value <= first(baseline) ? refinement.parameters : best.parameters
    loss, residue = pp_training_objective(pieces, parameters, shape, protocol; refined=true)
    boundary =
        any((parameters .- lower) ./ width .< 1e-4) ||
        any((upper .- parameters) ./ width .< 1e-4)
    return (; parameters, residue, loss, boundary, results, refinement, shape)
end

function fit_energy_shape(training, protocol, base)
    grid = protocol["energy_shape_grid"]
    settings = AdvancedVyazovkinConfig(
        protocol["energy_shape_delta_alpha"],
        base.integration_points,
        base.energy_bounds_kJ_per_mol,
        base.optimizer_grid_points,
        base.optimizer_tolerance_kJ_per_mol,
        base.optimizer_max_iterations,
        base.boundary_tolerance_kJ_per_mol,
        base.flat_relative_tolerance,
    )
    curves = TGConversionResult[observation.conversion for observation in training]
    result = tg_activation_energy(
        :advanced_vyazovkin,
        curves,
        grid;
        advanced=settings,
        time_subdivisions=protocol["energy_shape_subdivisions"],
    )
    all(p -> p.status === :estimated, result.points) ||
        error("training-only E(α) profile is unavailable for $(result.run_ids)")
    energy = Float64[p.energy_kJ_per_mol for p in result.points]
    reference_index = only(findall(==(0.5), grid))
    return (; shape=(grid, energy .- energy[reference_index]), energy)
end

function predict_pp(fit, program, protocol; refined=false)
    survival = reactive_history(
        program, program.history.time_min, fit.parameters, fit.shape, protocol; refined
    )
    return program.first_mass_mg .* (fit.residue .+ (1 - fit.residue) .* survival)
end

function score_pp(prediction, held, rows)
    errors = prediction[rows] .- held.mass_mg[rows]
    time = view(held.program.history.time_min, rows)
    w = weights_for_time(time)
    rmse = sqrt(sum(w .* abs2.(errors)))
    return Dict(
        "rows" => length(rows),
        "rmse_mg" => rmse,
        "rmse_percent_first_mass" => 100 * rmse / held.program.first_mass_mg,
        "bias_mg" => sum(w .* errors),
        "mae_mg" => sum(w .* abs.(errors)),
    )
end

function fold_record(id, held, training, model, fitted, shape_energy, protocol)
    train_ids = getfield.(training, :id)
    held.id in train_ids && error("held-out mass leaked into fitting")
    predicted = predict_pp(fitted, held.program, protocol)
    refined = predict_pp(fitted, held.program, protocol; refined=true)
    shift = maximum(abs, predicted .- refined)
    shift <= protocol["max_prediction_refinement_mg"] ||
        error("held-out integration refinement failed")
    return Dict{String,Any}(
        "fold" => id,
        "held_out_run" => held.id,
        "model" => String(model),
        "training_run_ids" => train_ids,
        "first_mass_mg" => held.program.first_mass_mg,
        "energy_at_half_conversion_kJ_per_mol" => fitted.parameters[1],
        "log_reference_rate_min_inv" => fitted.parameters[2],
        "reaction_order" => fitted.parameters[3],
        "retained_fraction" => fitted.residue,
        "training_loss" => fitted.loss,
        "parameter_boundary" => fitted.boundary,
        "training_refinement" =>
            Dict(string(k) => v for (k, v) in pairs(fitted.refinement)),
        "training_energy_profile_kJ_per_mol" => shape_energy,
        "starts" => [
            Dict(
                "parameters" => r.parameters,
                "loss" => r.value,
                "converged" => r.converged,
                "iterations" => r.iterations,
            ) for r in fitted.results
        ],
        "observed_final_mass_mg" => last(held.mass_mg),
        "predicted_final_mass_mg" => last(refined),
        "maximum_refinement_shift_mg" => shift,
        "full" => score_pp(refined, held, 1:length(refined)),
        "reaction_window" => score_pp(refined, held, held.reaction_rows),
        "predicted_mass_mg" => refined,
    )
end

function read_pp_observations(manifest, conversion_protocol, protocol)
    specs = filter(spec -> spec["group"] == "experimental_tg", manifest["runs"])
    return map(specs) do spec
        run = import_tg(spec, manifest)
        selected, _, reference = tg_case_curves(run, spec, conversion_protocol)
        time = time_min(selected)
        temperature = temperature_K(selected)
        mass = signal_values(selected, :mass)
        first(mass) > 0 || error("first measured mass is not positive")
        lo, hi = protocol["reaction_scoring_window_C"] .+ 273.15
        indices = findall(T -> lo <= T <= hi, temperature)
        rows = first(indices):last(indices)
        indices == collect(rows) || error("reaction scoring window is not contiguous")
        program = PPProgram(
            SampledTemperatureProgram(time .- first(time), temperature), first(mass)
        )
        return PPObservation(run.id, run.provenance.sha256, program, mass, rows, reference)
    end
end

function check_analytic_prediction(protocol)
    E, log_rate, n = (128.0, -1.3, 1.4)
    residue = 0.027
    reference = protocol["reference_temperature_K"]
    runs = map([595.0, 635.0, 695.0, 725.0]) do T
        k = exp(log_rate - E * (inv(T) - inv(reference)) / PP_R)
        t = collect(range(0.0, 18 / k; length=201))
        mass = [
            11 * (residue + (1 - residue) * (1 + (n - 1) * k * x)^(-1 / (n - 1))) for x in t
        ]
        program = PPProgram(SampledTemperatureProgram(t, fill(T, length(t))), first(mass))
        return PPObservation("analytic_$T", "", program, mass, 1:length(t), nothing)
    end
    fit = fit_pp(runs[1:3], nothing, protocol)
    prediction = predict_pp(fit, runs[4].program, protocol; refined=true)
    max_error = maximum(abs, prediction .- runs[4].mass_mg)
    max_error < 0.001 || error("analytic mass prediction control failed")
    return Dict(
        "truth" => [E, log_rate, n, residue],
        "fitted" => [fit.parameters; fit.residue],
        "held_out_max_error_mg" => max_error,
        "training_loss" => fit.loss,
    )
end

function audit_public_pp_prediction(output)
    protocol = TOML.parsefile(joinpath(TG_ROOT, PP_PROTOCOL))
    manifest = TOML.parsefile(joinpath(TG_ROOT, TG_IMPORT_MANIFEST))
    conversion = TOML.parsefile(joinpath(TG_ROOT, TG_PROTOCOL))
    base =
        load_config(joinpath(TG_ROOT, "config", "analysis_defaults.toml")).isoconversional.advanced_vyazovkin
    control = check_analytic_prediction(protocol)
    println("Analytic held-out control passed at $(control["held_out_max_error_mg"]) mg.")
    observations = read_pp_observations(manifest, conversion, protocol)
    omitted = protocol["omitted_sensitivity_run"]
    length(observations) == 5 && count(r -> r.id == omitted, observations) == 1 ||
        error("unexpected PP roster")
    folds = [
        ("all_rates_loocv", held, filter(r -> r.id != held.id, observations)) for
        held in observations
    ]
    append!(
        folds,
        [
            (
                "omit_10_loocv",
                held,
                filter(r -> r.id != held.id && r.id != omitted, observations),
            ) for held in observations if held.id != omitted
        ],
    )
    records = Vector{Vector{Dict{String,Any}}}(undef, length(folds))
    Threads.@threads for index in eachindex(folds)
        id, held, training = folds[index]
        profile = fit_energy_shape(training, protocol, base)
        records[index] = map(collect(PP_MODELS)) do model
            shape = model === :constant_energy ? nothing : profile.shape
            fit = fit_pp(training, shape, protocol)
            record = fold_record(
                id,
                held,
                training,
                model,
                fit,
                model === :constant_energy ? Float64[] : profile.energy,
                protocol,
            )
            println(
                "$id $(held.id) $model reaction RMSE $(round(record["reaction_window"]["rmse_mg"]; digits=5)) mg",
            )
            return record
        end
    end
    inputs = (PP_PROTOCOL, TG_IMPORT_MANIFEST, TG_PROTOCOL, "config/analysis_defaults.toml")
    report = Dict{String,Any}(
        "schema_version" => 1,
        "protocol" => PP_PROTOCOL,
        "control" => control,
        "input_sha256" => Dict(
            path => bytes2hex(open(sha256, joinpath(TG_ROOT, path))) for path in inputs
        ),
        "runs" => [
            Dict(
                "id" => r.id,
                "source_sha256" => r.sha256,
                "reaction_rows" => [first(r.reaction_rows), last(r.reaction_rows)],
            ) for r in observations
        ],
        "predictions" => reduce(vcat, records),
    )
    open(output, "w") do io
        return TOML.print(io, report; sorted=true)
    end
    println("Wrote $output")
    return report
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) <= 1 || error("Expected at most one report path")
    audit_public_pp_prediction(
        isempty(ARGS) ? joinpath(TG_ROOT, PP_REPORT) : abspath(only(ARGS))
    )
end
