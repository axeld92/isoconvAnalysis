"""
Raised when an M6 reaction-model or compensation calculation violates its scientific or
numerical contract.
"""
struct KineticModelError <: Exception
    context::String
    message::String
end

function Base.showerror(io::IO, error::KineticModelError)
    return print(io, "kinetic-model failure for '$(error.context)': $(error.message)")
end

"""
Named differential reaction model. `family` is descriptive and does not establish a physical
mechanism for an empirical fit.
"""
struct ReactionModelSpec
    name::Symbol
    label::String
    family::Symbol
    exponent_names::Vector{Symbol}
end

const _REACTION_MODEL_REGISTRY = ReactionModelSpec[
    ReactionModelSpec(:f1, "first order", :reaction_order, Symbol[]),
    ReactionModelSpec(:f2, "second order", :reaction_order, Symbol[]),
    ReactionModelSpec(:f3, "third order", :reaction_order, Symbol[]),
    ReactionModelSpec(:f4, "fourth order", :reaction_order, Symbol[]),
    ReactionModelSpec(:a2, "Avrami–Erofeev A2", :nucleation_growth, Symbol[]),
    ReactionModelSpec(:a3, "Avrami–Erofeev A3", :nucleation_growth, Symbol[]),
    ReactionModelSpec(:a4, "Avrami–Erofeev A4", :nucleation_growth, Symbol[]),
    ReactionModelSpec(:r2, "contracting cylinder R2", :geometrical_contraction, Symbol[]),
    ReactionModelSpec(:r3, "contracting sphere R3", :geometrical_contraction, Symbol[]),
    ReactionModelSpec(:d1, "one-dimensional diffusion D1", :diffusion, Symbol[]),
    ReactionModelSpec(:d2, "two-dimensional diffusion D2", :diffusion, Symbol[]),
    ReactionModelSpec(:d3, "Jander diffusion D3", :diffusion, Symbol[]),
    ReactionModelSpec(:d4, "Ginstling–Brounshtein diffusion D4", :diffusion, Symbol[]),
    ReactionModelSpec(:p2, "power law P2", :power_law, Symbol[]),
    ReactionModelSpec(:p3, "power law P3", :power_law, Symbol[]),
    ReactionModelSpec(:p4, "power law P4", :power_law, Symbol[]),
    ReactionModelSpec(:random_scission, "random scission", :random_scission, Symbol[]),
    ReactionModelSpec(
        :sestak_berggren_2, "truncated Šesták–Berggren", :empirical, [:m, :n]
    ),
    ReactionModelSpec(
        :sestak_berggren_3, "generalized Šesták–Berggren", :empirical, [:m, :n, :p]
    ),
]

"""
    reaction_model_registry()

Return the complete ordered M6 differential-model registry. Labels are conventional model
names, not mechanistic assignments for this dataset.
"""
function reaction_model_registry()
    return [
        ReactionModelSpec(spec.name, spec.label, spec.family, copy(spec.exponent_names)) for
        spec in _REACTION_MODEL_REGISTRY
    ]
end

function _reaction_spec(model::Symbol)
    index = findfirst(spec -> spec.name == model, _REACTION_MODEL_REGISTRY)
    isnothing(index) && throw(ArgumentError("unsupported reaction model $model"))
    return _REACTION_MODEL_REGISTRY[index]
end

function _validate_conversion(alpha::Real)
    value = Float64(alpha)
    isfinite(value) && 0 < value < 1 ||
        throw(DomainError(alpha, "reaction functions require 0 < alpha < 1"))
    return value
end

"""
    reaction_function(model, alpha; m=0, n=0, p=0)

Evaluate a dimensionless differential reaction function on `0 < alpha < 1`. The generalized
Šesták–Berggren convention is `alpha^m * (1-alpha)^n * [-log(1-alpha)]^p`;
the truncated convention fixes `p = 0`.
"""
function reaction_function(
    model::Symbol, alpha::Real; m::Real=0.0, n::Real=0.0, p::Real=0.0
)
    value = _validate_conversion(alpha)
    one_minus = 1 - value
    negative_log = -log1p(-value)
    result = if model == :f1
        one_minus
    elseif model == :f2
        one_minus^2
    elseif model == :f3
        one_minus^3
    elseif model == :f4
        one_minus^4
    elseif model == :a2
        2 * one_minus * negative_log^(1 / 2)
    elseif model == :a3
        3 * one_minus * negative_log^(2 / 3)
    elseif model == :a4
        4 * one_minus * negative_log^(3 / 4)
    elseif model == :r2
        2 * one_minus^(1 / 2)
    elseif model == :r3
        3 * one_minus^(2 / 3)
    elseif model == :d1
        inv(2 * value)
    elseif model == :d2
        inv(negative_log)
    elseif model == :d3
        (3 / 2) * one_minus^(2 / 3) / (1 - one_minus^(1 / 3))
    elseif model == :d4
        (3 / 2) / (one_minus^(-1 / 3) - 1)
    elseif model == :p2
        2 * value^(1 / 2)
    elseif model == :p3
        3 * value^(2 / 3)
    elseif model == :p4
        4 * value^(3 / 4)
    elseif model == :random_scission
        2 * (sqrt(value) - value)
    elseif model == :sestak_berggren_2
        value^Float64(m) * one_minus^Float64(n)
    elseif model == :sestak_berggren_3
        value^Float64(m) * one_minus^Float64(n) * negative_log^Float64(p)
    else
        throw(ArgumentError("unsupported reaction model $model"))
    end
    isfinite(result) && result > 0 || throw(
        DomainError(
            (model, alpha),
            "reaction function is not finite and positive at this conversion",
        ),
    )
    return Float64(result)
end

function reaction_function(
    model::Symbol, alpha::AbstractVector{<:Real}; m::Real=0.0, n::Real=0.0, p::Real=0.0
)
    return [reaction_function(model, value; m=m, n=n, p=p) for value in alpha]
end

"""
One coefficient estimate and its local two-sided Student-t interval. The interval is
conditional on the selected model, conversion interval, and independent-error approximation.
"""
struct KineticParameterEstimate
    name::Symbol
    unit::String
    estimate::Float64
    standard_error::Float64
    confidence_lower::Float64
    confidence_upper::Float64
end

"""
Pointwise validation evidence for one heating-rate experiment omitted from fitting.
"""
struct HeatingRateValidation
    experiment_id::String
    point_count::Int
    log_rate_rmse::Float64
    relative_rate_rmse::Float64
    mean_log_rate_bias::Float64
end

"""
Fit, information-criterion, identifiability, and leave-one-heating-rate-out evidence for one
kinetic-triplet candidate.
"""
struct KineticTripletDiagnostics
    status::Symbol
    point_count::Int
    experiment_count::Int
    excluded_nonpositive_rate_count::Int
    degrees_of_freedom::Int
    design_rank::Int
    parameter_count::Int
    normalized_design_condition::Float64
    maximum_absolute_parameter_correlation::Float64
    log_residual_sum_squares::Float64
    log_rate_rmse::Float64
    log_rate_r_squared::Float64
    relative_rate_rmse::Float64
    median_durbin_watson::Float64
    aicc::Float64
    bic::Float64
    cross_validation_log_rate_rmse::Float64
    cross_validation_relative_rate_rmse::Float64
    warnings::Vector{String}
end

"""
One overall-conversion kinetic triplet fitted jointly across heating programs. The
pre-exponential factor uses `min^-1` because the package time boundary is minutes.
"""
struct KineticTripletResult
    model::ReactionModelSpec
    composition::Composition
    experiment_ids::Vector{String}
    point_experiment_ids::Vector{String}
    alpha::Vector{Float64}
    temperature_K::Vector{Float64}
    observed_rate_per_min::Vector{Float64}
    fitted_rate_per_min::Vector{Float64}
    log_rate_residual::Vector{Float64}
    parameters::Vector{KineticParameterEstimate}
    parameter_covariance::Matrix{Float64}
    parameter_correlation::Matrix{Float64}
    validations::Vector{HeatingRateValidation}
    configuration::ReactionModelConfig
    diagnostics::KineticTripletDiagnostics
    analysis_fingerprint::String
end

"""
Complete M6 comparison. Every candidate, the raw criterion minimum, structural eligibility,
and the final selection remain visible.
"""
struct KineticModelComparison
    composition::Composition
    experiment_ids::Vector{String}
    criterion::Symbol
    criterion_minimum_model::Symbol
    selected_model::Symbol
    criterion_values::Vector{Float64}
    criterion_deltas::Vector{Float64}
    structurally_eligible::Vector{Bool}
    results::Vector{KineticTripletResult}
    status::Symbol
    warnings::Vector{String}
    analysis_fingerprint::String
end

"""
Single-heating-rate invariant-parameter/compensation diagnostic. A high compensation-line R²
does not establish a physical compensation effect or validate single-rate kinetic parameters.
"""
struct CompensationResult
    experiment_id::String
    model_names::Vector{Symbol}
    model_activation_energies_kJ_per_mol::Vector{Float64}
    model_log_preexponentials_per_min::Vector{Float64}
    reference_activation_energy_kJ_per_mol::Float64
    estimated_log_preexponential_per_min::Float64
    estimated_preexponential_per_min::Float64
    confidence_lower_preexponential_per_min::Float64
    confidence_upper_preexponential_per_min::Float64
    regression::LinearRegressionDiagnostics
    warnings::Vector{String}
    configuration::ReactionModelConfig
    analysis_fingerprint::String
end

struct _ReactionFitData
    experiment_ids::Vector{String}
    point_experiment_ids::Vector{String}
    alpha::Vector{Float64}
    temperature_K::Vector{Float64}
    rate_per_min::Vector{Float64}
    excluded_nonpositive_rate_count::Int
end

function _kinetic_error(context, message)
    return throw(KineticModelError(String(context), String(message)))
end

function _selected_reaction_points(
    experiment::ProcessedExperiment, config::ReactionModelConfig
)
    lower_alpha, upper_alpha = config.fit_conversion_range
    lower_crossings = _crossing_indices(experiment.analysis_alpha, lower_alpha)
    upper_crossings = _crossing_indices(experiment.analysis_alpha, upper_alpha)
    isempty(lower_crossings) && _kinetic_error(
        experiment.source_id, "lower conversion $lower_alpha has no upward crossing"
    )
    isempty(upper_crossings) && _kinetic_error(
        experiment.source_id, "upper conversion $upper_alpha has no upward crossing"
    )
    lower_index = first(lower_crossings)
    upper_candidates = filter(index -> index >= lower_index, upper_crossings)
    isempty(upper_candidates) && _kinetic_error(
        experiment.source_id, "upper conversion crossing precedes the lower crossing"
    )
    upper_index = first(upper_candidates)
    lower_sample = interpolate_at_conversion(experiment, lower_alpha)
    upper_sample = interpolate_at_conversion(experiment, upper_alpha)
    interior = (lower_index + 1):upper_index
    alpha = vcat(lower_alpha, experiment.analysis_alpha[interior], upper_alpha)
    temperature = vcat(
        lower_sample.temperature_K,
        experiment.temperature_K[interior],
        upper_sample.temperature_K,
    )
    rate = vcat(
        lower_sample.dalpha_dt_min_inv,
        experiment.dalpha_dt_min_inv[interior],
        upper_sample.dalpha_dt_min_inv,
    )
    all(>(0), diff(temperature)) || _kinetic_error(
        experiment.source_id, "selected temperatures are not strictly increasing"
    )
    valid = [
        isfinite(alpha[index]) &&
            isfinite(temperature[index]) &&
            isfinite(rate[index]) &&
            rate[index] > config.minimum_rate_per_min for index in eachindex(alpha)
    ]
    excluded = count(!, valid)
    valid_indices = findall(valid)
    length(valid_indices) >= 10 || _kinetic_error(
        experiment.source_id,
        "fewer than ten finite rates exceed $(config.minimum_rate_per_min) min^-1",
    )
    retained = valid_indices[_fit_indices(
        length(valid_indices), config.maximum_points_per_experiment
    )]
    return (
        alpha=alpha[retained],
        temperature_K=temperature[retained],
        rate_per_min=rate[retained],
        excluded=excluded,
    )
end

function _prepare_reaction_data(experiments, config)
    ordered = sort(collect(experiments); by=experiment -> experiment.source_id)
    curves = [_selected_reaction_points(experiment, config) for experiment in ordered]
    return _ReactionFitData(
        [experiment.source_id for experiment in ordered],
        reduce(
            vcat,
            [
                fill(experiment.source_id, length(curve.alpha)) for
                (experiment, curve) in zip(ordered, curves)
            ],
        ),
        reduce(vcat, getproperty.(curves, :alpha)),
        reduce(vcat, getproperty.(curves, :temperature_K)),
        reduce(vcat, getproperty.(curves, :rate_per_min)),
        sum(getproperty.(curves, :excluded)),
    )
end

function _validate_reaction_group(experiments, config)
    length(experiments) >= config.minimum_experiments || _kinetic_error(
        "group",
        "at least $(config.minimum_experiments) experiments are required; received $(length(experiments))",
    )
    length(unique(experiment.source_id for experiment in experiments)) ==
    length(experiments) || _kinetic_error("group", "experiment identifiers must be unique")
    first_composition = first(experiments).composition
    all(experiment -> experiment.composition == first_composition, experiments) ||
        _kinetic_error("group", "all experiments must have the same composition")
    first_initial = first(experiments).initial_temperature_K
    first_final = first(experiments).final_temperature_K
    all(
        experiment ->
            experiment.initial_temperature_K == first_initial &&
            experiment.final_temperature_K == first_final,
        experiments,
    ) || _kinetic_error(
        "group", "all experiments must use the same conversion-reference temperatures"
    )
    return nothing
end

function _design(model::Symbol, data::_ReactionFitData)
    inverse_arrhenius = @. -1 / (GAS_CONSTANT_KJ_PER_MOL_K * data.temperature_K)
    if model in _FIXED_REACTION_MODELS
        function_values = reaction_function(model, data.alpha)
        response = log.(data.rate_per_min) .- log.(function_values)
        return hcat(ones(length(response)), inverse_arrhenius),
        response, log.(function_values),
        [:log_preexponential_per_min, :activation_energy_kJ_per_mol]
    elseif model == :sestak_berggren_2
        response = log.(data.rate_per_min)
        design = hcat(
            ones(length(response)), inverse_arrhenius, log.(data.alpha), log1p.(-data.alpha)
        )
        return design,
        response, zeros(length(response)),
        [:log_preexponential_per_min, :activation_energy_kJ_per_mol, :m, :n]
    elseif model == :sestak_berggren_3
        response = log.(data.rate_per_min)
        design = hcat(
            ones(length(response)),
            inverse_arrhenius,
            log.(data.alpha),
            log1p.(-data.alpha),
            log.(-log1p.(-data.alpha)),
        )
        return design,
        response, zeros(length(response)),
        [:log_preexponential_per_min, :activation_energy_kJ_per_mol, :m, :n, :p]
    end
    return _kinetic_error(model, "model is not supported")
end

function _matrix_correlation(covariance)
    scales = sqrt.(max.(diag(covariance), 0.0))
    correlation = zeros(Float64, size(covariance))
    for row in axes(covariance, 1), column in axes(covariance, 2)
        denominator = scales[row] * scales[column]
        correlation[row, column] =
            denominator > 0 ? covariance[row, column] / denominator : 0.0
    end
    for index in axes(correlation, 1)
        correlation[index, index] = 1.0
    end
    return correlation
end

function _normalized_condition(design)
    normalized = copy(design)
    for column in axes(normalized, 2)
        column_norm = norm(view(normalized, :, column))
        column_norm > 0 && (normalized[:, column] ./= column_norm)
    end
    singular_values = svdvals(normalized)
    return if isempty(singular_values) || last(singular_values) <= 0
        Inf
    else
        first(singular_values) / last(singular_values)
    end
end

function _maximum_off_diagonal(correlation)
    values = [
        abs(correlation[row, column]) for
        row in axes(correlation, 1), column in axes(correlation, 2) if row != column
    ]
    return isempty(values) ? 0.0 : maximum(values)
end

function _parameter_unit(name::Symbol)
    name == :log_preexponential_per_min && return "ln(min^-1)"
    name == :activation_energy_kJ_per_mol && return "kJ/mol"
    return "1"
end

function _core_fit(model::Symbol, data::_ReactionFitData, config::ReactionModelConfig)
    design, response, offset, parameter_names = _design(model, data)
    point_count, parameter_count = size(design)
    point_count > parameter_count || _kinetic_error(model, "too few points for model")
    coefficients = design \ response
    fitted_log_rate = offset .+ design * coefficients
    observed_log_rate = log.(data.rate_per_min)
    residuals = observed_log_rate .- fitted_log_rate
    residual_sum_squares = max(sum(abs2, residuals), 0.0)
    degrees_of_freedom = point_count - parameter_count
    residual_variance = residual_sum_squares / degrees_of_freedom
    covariance =
        residual_variance .* pinv(transpose(design) * design; rtol=sqrt(eps(Float64)))
    standard_errors = sqrt.(max.(diag(covariance), 0.0))
    critical = quantile(TDist(degrees_of_freedom), 0.5 + config.confidence_level / 2)
    lower = coefficients .- critical .* standard_errors
    upper = coefficients .+ critical .* standard_errors
    parameter_estimates = [
        KineticParameterEstimate(
            parameter_names[index],
            _parameter_unit(parameter_names[index]),
            coefficients[index],
            standard_errors[index],
            lower[index],
            upper[index],
        ) for index in eachindex(parameter_names)
    ]
    correlation = _matrix_correlation(covariance)
    fitted_rate = exp.(fitted_log_rate)
    rate_residuals = data.rate_per_min .- fitted_rate
    total_log_sum_squares = sum(abs2, observed_log_rate .- mean(observed_log_rate))
    log_r_squared =
        total_log_sum_squares > 0 ? 1 - residual_sum_squares / total_log_sum_squares : NaN
    relative_rate_rmse = sqrt(
        sum(abs2, rate_residuals) / max(sum(abs2, data.rate_per_min), eps(Float64))
    )
    durbin_watson_values = Float64[]
    for experiment_id in data.experiment_ids
        indices = findall(==(experiment_id), data.point_experiment_ids)
        curve_residuals = residuals[indices]
        denominator = sum(abs2, curve_residuals)
        push!(
            durbin_watson_values,
            denominator > 0 ? sum(abs2, diff(curve_residuals)) / denominator : 2.0,
        )
    end
    aic =
        point_count * log(max(residual_sum_squares / point_count, eps(Float64))) +
        2 * parameter_count
    aicc = if point_count > parameter_count + 1
        aic + 2 * parameter_count * (parameter_count + 1) / (point_count - parameter_count - 1)
    else
        Inf
    end
    bic =
        point_count * log(max(residual_sum_squares / point_count, eps(Float64))) +
        parameter_count * log(point_count)
    return (
        coefficients=coefficients,
        parameter_names=parameter_names,
        parameter_estimates=parameter_estimates,
        covariance=covariance,
        correlation=correlation,
        fitted_log_rate=fitted_log_rate,
        fitted_rate=fitted_rate,
        residuals=residuals,
        degrees_of_freedom=degrees_of_freedom,
        rank=rank(design; rtol=sqrt(eps(Float64))),
        condition=_normalized_condition(design),
        maximum_correlation=_maximum_off_diagonal(correlation),
        rss=residual_sum_squares,
        log_rmse=sqrt(residual_sum_squares / point_count),
        log_r_squared=log_r_squared,
        relative_rate_rmse=relative_rate_rmse,
        median_durbin_watson=median(durbin_watson_values),
        aicc=aicc,
        bic=bic,
    )
end

function _parameter_value(core, name::Symbol)
    index = findfirst(==(name), core.parameter_names)
    isnothing(index) && return 0.0
    return core.coefficients[index]
end

function _predict_log_rate(core, model::Symbol, temperature_K, alpha)
    temperatures = Float64.(temperature_K)
    conversions = Float64.(alpha)
    length(temperatures) == length(conversions) ||
        throw(DimensionMismatch("temperature and conversion lengths differ"))
    energy = _parameter_value(core, :activation_energy_kJ_per_mol)
    log_preexponential = _parameter_value(core, :log_preexponential_per_min)
    m = _parameter_value(core, :m)
    n = _parameter_value(core, :n)
    p = _parameter_value(core, :p)
    log_function = log.(reaction_function(model, conversions; m=m, n=n, p=p))
    return @. log_preexponential - energy / (GAS_CONSTANT_KJ_PER_MOL_K * temperatures) +
        log_function
end

function _cross_validate(model, experiments, config)
    ordered = sort(collect(experiments); by=experiment -> experiment.source_id)
    validations = HeatingRateValidation[]
    all_residuals = Float64[]
    squared_rate_errors = 0.0
    squared_observed_rates = 0.0
    for held_out_index in eachindex(ordered)
        training = [
            experiment for
            (index, experiment) in enumerate(ordered) if index != held_out_index
        ]
        training_data = _prepare_reaction_data(training, config)
        core = _core_fit(model, training_data, config)
        held_out = ordered[held_out_index]
        test_data = _prepare_reaction_data([held_out], config)
        predicted_log = _predict_log_rate(
            core, model, test_data.temperature_K, test_data.alpha
        )
        observed_log = log.(test_data.rate_per_min)
        residuals = observed_log .- predicted_log
        predicted_rate = exp.(predicted_log)
        rate_errors = test_data.rate_per_min .- predicted_rate
        append!(all_residuals, residuals)
        squared_rate_errors += sum(abs2, rate_errors)
        squared_observed_rates += sum(abs2, test_data.rate_per_min)
        push!(
            validations,
            HeatingRateValidation(
                held_out.source_id,
                length(residuals),
                sqrt(mean(abs2, residuals)),
                sqrt(
                    sum(abs2, rate_errors) /
                    max(sum(abs2, test_data.rate_per_min), eps(Float64)),
                ),
                mean(residuals),
            ),
        )
    end
    return validations,
    sqrt(mean(abs2, all_residuals)),
    sqrt(squared_rate_errors / max(squared_observed_rates, eps(Float64)))
end

function _reaction_fingerprint(experiments, model, config)
    ordered = sort(collect(experiments); by=experiment -> experiment.source_id)
    fields = (
        model,
        [
            (experiment.source_id, experiment.source_sha256, experiment.config_fingerprint)
            for experiment in ordered
        ],
        config.models,
        config.compensation_models,
        config.fit_conversion_range,
        config.maximum_points_per_experiment,
        config.minimum_rate_per_min,
        config.selection_criterion,
        config.energy_bounds_kJ_per_mol,
        config.sestak_berggren_exponent_bounds,
        config.minimum_experiments,
        config.confidence_level,
        config.identifiability_condition_warning,
        config.parameter_correlation_warning,
        config.ambiguity_criterion_delta,
        config.maximum_cross_validation_log_rmse,
    )
    return bytes2hex(sha256(join(repr.(fields), "|")))
end

function _structurally_eligible(core, model, config)
    energy = _parameter_value(core, :activation_energy_kJ_per_mol)
    energy_valid =
        config.energy_bounds_kJ_per_mol[1] <= energy <= config.energy_bounds_kJ_per_mol[2]
    exponent_valid = all(
        name -> begin
            value = _parameter_value(core, name)
            config.sestak_berggren_exponent_bounds[1] <=
            value <=
            config.sestak_berggren_exponent_bounds[2]
        end,
        _reaction_spec(model).exponent_names,
    )
    return core.rank == length(core.parameter_names) &&
           core.condition <= config.identifiability_condition_warning &&
           energy_valid &&
           exponent_valid
end

function _build_kinetic_result(experiments, model, config)
    data = _prepare_reaction_data(experiments, config)
    core = _core_fit(model, data, config)
    validations, cv_log_rmse, cv_rate_rmse = _cross_validate(model, experiments, config)
    warnings = String[]
    data.excluded_nonpositive_rate_count > 0 && push!(
        warnings,
        "$(data.excluded_nonpositive_rate_count) nonpositive or nonfinite rate points were excluded because logarithms are undefined",
    )
    core.rank < length(core.parameter_names) &&
        push!(warnings, "design matrix is rank deficient")
    core.condition > config.identifiability_condition_warning && push!(
        warnings,
        "normalized design condition $(core.condition) exceeds the configured threshold",
    )
    core.maximum_correlation > config.parameter_correlation_warning && push!(
        warnings,
        "maximum absolute parameter correlation $(core.maximum_correlation) exceeds the configured threshold",
    )
    energy = _parameter_value(core, :activation_energy_kJ_per_mol)
    !(config.energy_bounds_kJ_per_mol[1] <= energy <= config.energy_bounds_kJ_per_mol[2]) &&
        push!(warnings, "activation energy lies outside the configured physical range")
    for exponent in _reaction_spec(model).exponent_names
        value = _parameter_value(core, exponent)
        !(
            config.sestak_berggren_exponent_bounds[1] <=
            value <=
            config.sestak_berggren_exponent_bounds[2]
        ) && push!(warnings, "$exponent lies outside the configured empirical range")
    end
    cv_log_rmse > config.maximum_cross_validation_log_rmse && push!(
        warnings,
        "leave-one-heating-rate-out log-rate RMSE exceeds the configured predictive threshold",
    )
    core.median_durbin_watson < 0.5 && push!(
        warnings,
        "strong positive residual autocorrelation makes pointwise intervals and information criteria conditional",
    )
    eligible = _structurally_eligible(core, model, config)
    status =
        if core.rank < length(core.parameter_names) ||
            core.condition > config.identifiability_condition_warning
            :nonidentifiable
        elseif !eligible
            :nonphysical
        elseif isempty(warnings)
            :ok
        else
            :warning
        end
    diagnostics = KineticTripletDiagnostics(
        status,
        length(data.alpha),
        length(data.experiment_ids),
        data.excluded_nonpositive_rate_count,
        core.degrees_of_freedom,
        core.rank,
        length(core.parameter_names),
        core.condition,
        core.maximum_correlation,
        core.rss,
        core.log_rmse,
        core.log_r_squared,
        core.relative_rate_rmse,
        core.median_durbin_watson,
        core.aicc,
        core.bic,
        cv_log_rmse,
        cv_rate_rmse,
        warnings,
    )
    ordered = sort(collect(experiments); by=experiment -> experiment.source_id)
    return KineticTripletResult(
        _reaction_spec(model),
        first(ordered).composition,
        data.experiment_ids,
        data.point_experiment_ids,
        data.alpha,
        data.temperature_K,
        data.rate_per_min,
        core.fitted_rate,
        core.residuals,
        core.parameter_estimates,
        core.covariance,
        core.correlation,
        validations,
        config,
        diagnostics,
        _reaction_fingerprint(ordered, model, config),
    )
end

"""
    fit_kinetic_triplet(experiments, model; configuration)

Fit one overall-conversion Arrhenius/reaction-model triplet jointly across heating programs.
The regression is performed in natural-log rate space with a single `ln(A)` intercept. Each
experiment contributes at most the configured number of deterministic points, and
leave-one-heating-rate-out validation is returned with the fit.
"""
function fit_kinetic_triplet(
    experiments::AbstractVector{<:ProcessedExperiment},
    model::Symbol;
    configuration::ReactionModelConfig,
)
    model in configuration.models || _kinetic_error(
        model, "model is not configured; allowed values are $(configuration.models)"
    )
    _validate_reaction_group(experiments, configuration)
    return _build_kinetic_result(experiments, model, configuration)
end

function fit_kinetic_triplet(
    experiments::AbstractVector{<:ProcessedExperiment},
    model::Symbol,
    config::AnalysisConfig,
)
    return fit_kinetic_triplet(experiments, model; configuration=config.reaction_models)
end

"""
    compare_reaction_models(experiments; configuration)

Fit every configured reaction model. Selection uses BIC or AICc among candidates with full
rank, acceptable condition, and configured energy/exponent ranges. Predictive failure and
criterion ambiguity remain explicit statuses rather than deleting the result.
"""
function compare_reaction_models(
    experiments::AbstractVector{<:ProcessedExperiment}; configuration::ReactionModelConfig
)
    _validate_reaction_group(experiments, configuration)
    results = [
        _build_kinetic_result(experiments, model, configuration) for
        model in configuration.models
    ]
    values = if configuration.selection_criterion == :bic
        getproperty.(getproperty.(results, :diagnostics), :bic)
    else
        getproperty.(getproperty.(results, :diagnostics), :aicc)
    end
    criterion_index = argmin(values)
    deltas = values .- values[criterion_index]
    eligible = [
        result.diagnostics.design_rank == result.diagnostics.parameter_count &&
            result.diagnostics.normalized_design_condition <=
            configuration.identifiability_condition_warning &&
            configuration.energy_bounds_kJ_per_mol[1] <=
            only(
                parameter.estimate for parameter in result.parameters if
                parameter.name == :activation_energy_kJ_per_mol
            ) <=
            configuration.energy_bounds_kJ_per_mol[2] &&
            all(
                parameter -> if parameter.name in (:m, :n, :p)
                    configuration.sestak_berggren_exponent_bounds[1] <=
                    parameter.estimate <=
                    configuration.sestak_berggren_exponent_bounds[2]
                else
                    true
                end,
                result.parameters,
            ) for result in results
    ]
    eligible_indices = findall(eligible)
    warnings = String[]
    selected_index = if isempty(eligible_indices)
        push!(warnings, "no candidate passed structural eligibility; raw minimum retained")
        criterion_index
    else
        eligible_indices[argmin(values[eligible_indices])]
    end
    constraint_filtered = selected_index != criterion_index
    ambiguity_values = isempty(eligible_indices) ? values : values[eligible_indices]
    ordered_ambiguity_deltas = sort(ambiguity_values .- minimum(ambiguity_values))
    criterion_ambiguous =
        length(ordered_ambiguity_deltas) > 1 &&
        ordered_ambiguity_deltas[2] < configuration.ambiguity_criterion_delta
    predictive_failure =
        results[selected_index].diagnostics.cross_validation_log_rate_rmse >
        configuration.maximum_cross_validation_log_rmse
    constraint_filtered &&
        push!(warnings, "raw criterion minimum failed structural eligibility")
    criterion_ambiguous &&
        push!(warnings, "best and second-best models are criterion-ambiguous")
    predictive_failure &&
        push!(warnings, "selected model fails the cross-validation error gate")
    status = if isempty(eligible_indices)
        :warning
    elseif predictive_failure
        :predictive_warning
    elseif constraint_filtered
        :constraint_filtered
    elseif criterion_ambiguous
        :ambiguous
    else
        :selected
    end
    ordered = sort(collect(experiments); by=experiment -> experiment.source_id)
    fingerprint = bytes2hex(
        sha256(
            join(
                [
                    _reaction_fingerprint(ordered, model, configuration) for
                    model in configuration.models
                ],
                "|",
            ),
        ),
    )
    return KineticModelComparison(
        first(ordered).composition,
        [experiment.source_id for experiment in ordered],
        configuration.selection_criterion,
        results[criterion_index].model.name,
        results[selected_index].model.name,
        values,
        deltas,
        eligible,
        results,
        status,
        warnings,
        fingerprint,
    )
end

function compare_reaction_models(
    experiments::AbstractVector{<:ProcessedExperiment}, config::AnalysisConfig
)
    return compare_reaction_models(experiments; configuration=config.reaction_models)
end

function _result_parameter(result::KineticTripletResult, name::Symbol)
    parameter = findfirst(value -> value.name == name, result.parameters)
    isnothing(parameter) && throw(ArgumentError("result has no parameter $name"))
    return result.parameters[parameter].estimate
end

function _optional_result_parameter(
    result::KineticTripletResult, name::Symbol, default::Float64=0.0
)
    parameter = findfirst(value -> value.name == name, result.parameters)
    return isnothing(parameter) ? default : result.parameters[parameter].estimate
end

"""
    predict_rate(result, temperature_K, alpha)

Evaluate the fitted central conversion rate in `min^-1`. This is a pointwise rate evaluation,
not an integrated conversion-history prediction; integration and state enforcement belong to
M7.
"""
function predict_rate(result::KineticTripletResult, temperature_K::Real, alpha::Real)
    temperature = Float64(temperature_K)
    isfinite(temperature) && temperature > 0 ||
        throw(DomainError(temperature_K, "temperature must be positive kelvin"))
    log_preexponential = _result_parameter(result, :log_preexponential_per_min)
    energy = _result_parameter(result, :activation_energy_kJ_per_mol)
    m = _optional_result_parameter(result, :m)
    n = _optional_result_parameter(result, :n)
    p = _optional_result_parameter(result, :p)
    function_value = reaction_function(result.model.name, alpha; m=m, n=n, p=p)
    return exp(log_preexponential - energy / (GAS_CONSTANT_KJ_PER_MOL_K * temperature)) *
           function_value
end

function _prediction_vector(result::KineticTripletResult, temperature_K, alpha)
    temperature = Float64(temperature_K)
    conversion = _validate_conversion(alpha)
    row = Float64[1.0, -1 / (GAS_CONSTANT_KJ_PER_MOL_K * temperature)]
    if result.model.name == :sestak_berggren_2
        append!(row, (log(conversion), log1p(-conversion)))
    elseif result.model.name == :sestak_berggren_3
        append!(row, (log(conversion), log1p(-conversion), log(-log1p(-conversion))))
    end
    return row
end

"""
    predict_rate_confidence_interval(result, temperature_K, alpha)

Propagate the local coefficient covariance to a two-sided confidence interval for the mean
log rate and transform it to `min^-1`. This is not a future-observation prediction interval.
"""
function predict_rate_confidence_interval(
    result::KineticTripletResult, temperature_K::Real, alpha::Real
)
    row = _prediction_vector(result, temperature_K, alpha)
    central = predict_rate(result, temperature_K, alpha)
    log_standard_error = sqrt(max(dot(row, result.parameter_covariance * row), 0.0))
    critical = quantile(
        TDist(result.diagnostics.degrees_of_freedom),
        0.5 + result.configuration.confidence_level / 2,
    )
    half_width = critical * log_standard_error
    return central * exp(-half_width), central, central * exp(half_width)
end

"""
    analyze_compensation(experiment, reference_energy; configuration)

Fit every configured fixed reaction model to one heating-rate experiment, regress its
`ln(A_i)` values against `E_i`, and evaluate the compensation line at a supplied independent
activation energy. The output is explicitly diagnostic because single-rate Arrhenius/model
separation is ambiguous and the parameter pairs share the same data.
"""
function analyze_compensation(
    experiment::ProcessedExperiment,
    reference_energy_kJ_per_mol::Real;
    configuration::ReactionModelConfig,
)
    reference_energy = Float64(reference_energy_kJ_per_mol)
    isfinite(reference_energy) && reference_energy > 0 || _kinetic_error(
        experiment.source_id, "reference activation energy must be finite and positive"
    )
    data = _prepare_reaction_data([experiment], configuration)
    cores = [
        _core_fit(model, data, configuration) for model in configuration.compensation_models
    ]
    energies = [_parameter_value(core, :activation_energy_kJ_per_mol) for core in cores]
    log_preexponentials = [
        _parameter_value(core, :log_preexponential_per_min) for core in cores
    ]
    regression = _linear_regression(energies, log_preexponentials)
    estimated_log = regression.intercept + regression.slope * reference_energy
    centered = energies .- mean(energies)
    prediction_standard_error =
        regression.residual_standard_error * sqrt(
            1 / length(energies) +
            (reference_energy - mean(energies))^2 / dot(centered, centered),
        )
    critical = quantile(
        TDist(regression.degrees_of_freedom), 0.5 + configuration.confidence_level / 2
    )
    lower_log = estimated_log - critical * prediction_standard_error
    upper_log = estimated_log + critical * prediction_standard_error
    warnings = [
        "compensation pairs come from the same single-rate data; high R² is not evidence of a physical compensation law",
    ]
    !(minimum(energies) <= reference_energy <= maximum(energies)) && push!(
        warnings,
        "reference activation energy lies outside the fitted compensation range",
    )
    any(
        energy -> !(
            configuration.energy_bounds_kJ_per_mol[1] <=
            energy <=
            configuration.energy_bounds_kJ_per_mol[2]
        ),
        energies,
    ) && push!(
        warnings, "one or more single-rate model energies are outside configured bounds"
    )
    fields = (
        experiment.source_id,
        experiment.source_sha256,
        experiment.config_fingerprint,
        reference_energy,
        configuration.compensation_models,
        configuration.fit_conversion_range,
        configuration.maximum_points_per_experiment,
        configuration.minimum_rate_per_min,
        configuration.confidence_level,
    )
    return CompensationResult(
        experiment.source_id,
        copy(configuration.compensation_models),
        energies,
        log_preexponentials,
        reference_energy,
        estimated_log,
        exp(estimated_log),
        exp(lower_log),
        exp(upper_log),
        regression,
        warnings,
        configuration,
        bytes2hex(sha256(join(repr.(fields), "|"))),
    )
end

function analyze_compensation(
    experiment::ProcessedExperiment,
    reference_energy_kJ_per_mol::Real,
    config::AnalysisConfig,
)
    return analyze_compensation(
        experiment, reference_energy_kJ_per_mol; configuration=config.reaction_models
    )
end
