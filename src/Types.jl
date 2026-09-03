"""
Dry-mass fractions of the two feedstocks in one sample.
"""
struct Composition
    waste_tire::Float64
    acrocomia_endocarp::Float64

    function Composition(waste_tire::Real, acrocomia_endocarp::Real)
        tire = Float64(waste_tire)
        endocarp = Float64(acrocomia_endocarp)
        all(isfinite, (tire, endocarp)) ||
            throw(ArgumentError("composition fractions must be finite"))
        all(fraction -> 0.0 <= fraction <= 1.0, (tire, endocarp)) ||
            throw(ArgumentError("composition fractions must lie in [0, 1]"))
        isapprox(tire + endocarp, 1.0; atol=1.0e-10, rtol=0.0) ||
            throw(ArgumentError("composition fractions must sum to one"))
        return new(tire, endocarp)
    end
end

Composition(waste_tire::Real) = Composition(waste_tire, 1.0 - waste_tire)

"""
One validation finding attached to a raw experiment.
"""
struct ValidationIssue
    severity::Symbol
    code::Symbol
    message::String
end

function Base.show(io::IO, issue::ValidationIssue)
    return print(io, "$(issue.severity):$(issue.code): $(issue.message)")
end

"""
Raised when an `Experiment` violates a structural data contract.
"""
struct ExperimentValidationError <: Exception
    experiment_id::String
    issues::Vector{ValidationIssue}
end

function Base.showerror(io::IO, error::ExperimentValidationError)
    print(io, "invalid experiment '$(error.experiment_id)'")
    for issue in error.issues
        print(io, "\n  - ")
        show(io, issue)
    end
end

"""
Raised when a dataset catalog is missing required or valid metadata.
"""
struct DatasetConfigError <: Exception
    path::String
    message::String
end

function Base.showerror(io::IO, error::DatasetConfigError)
    return print(io, "invalid dataset catalog at $(error.path): $(error.message)")
end

"""
Raised when an explicitly requested source variable cannot be imported safely.
"""
struct DataImportError <: Exception
    source_file::String
    source_variable::String
    message::String
end

function Base.showerror(io::IO, error::DataImportError)
    variable = isempty(error.source_variable) ? "" : " variable '$(error.source_variable)'"
    return print(
        io, "data import failed for $(error.source_file)$variable: $(error.message)"
    )
end

"""
Declarative metadata for one experiment within a MAT source.
"""
struct ExperimentSpec
    id::String
    source_variable::String
    sample_id::String
    composition::Composition
    nominal_heating_rate_K_per_min::Float64
    nominal_hold_temperature_celsius::Union{Nothing,Float64}
    initial_mass_mg::Union{Missing,Float64}
    acquisition_datetime::Union{Missing,DateTime}
    raw_export_file::Union{Nothing,String}
    legacy_alias::Union{Nothing,String}
    replicate::Union{Missing,String}
    notes::Vector{String}
end

"""
One MAT file and the experiments explicitly selected from it.
"""
struct DatasetSourceSpec
    id::String
    role::Symbol
    kind::Symbol
    format::Symbol
    source_file::String
    expected_sha256::String
    atmosphere::String
    required_columns::Vector{Symbol}
    optional_columns::Vector{Symbol}
    notes::Vector{String}
    experiments::Vector{ExperimentSpec}
end

"""
Validated, path-resolved catalog of selected experimental sources.
"""
struct DatasetCatalog
    source_path::String
    schema_version::Int
    sources::Vector{DatasetSourceSpec}
end

"""
    Experiment

Raw imported experiment in canonical boundary units. Rows are preserved exactly, including
recognized invalid acquisition rows. `import_warnings` records those conditions so M3 can
apply an explicit preprocessing policy.
"""
struct Experiment
    id::String
    role::Symbol
    kind::Symbol
    source_file::String
    source_variable::String
    source_sha256::String
    sample_id::String
    composition::Composition
    temperature_K::Vector{Float64}
    time_min::Vector{Float64}
    mass_percent::Vector{Float64}
    dsc_mW_per_mg::Union{Nothing,Vector{Float64}}
    purge_flow_mL_per_min::Union{Nothing,Vector{Float64}}
    protective_flow_mL_per_min::Union{Nothing,Vector{Float64}}
    sensitivity_uV_per_mW::Union{Nothing,Vector{Float64}}
    segment::Union{Nothing,Vector{Union{Missing,Int}}}
    nominal_heating_rate_K_per_min::Float64
    nominal_hold_temperature_celsius::Union{Nothing,Float64}
    atmosphere::String
    initial_mass_mg::Union{Missing,Float64}
    acquisition_datetime::Union{Missing,DateTime}
    raw_export_file::Union{Nothing,String}
    legacy_alias::Union{Nothing,String}
    replicate::Union{Missing,String}
    notes::Vector{String}
    import_warnings::Vector{String}

    function Experiment(
        ::Val{:validated},
        id,
        role,
        kind,
        source_file,
        source_variable,
        source_sha256,
        sample_id,
        composition,
        temperature_K,
        time_min,
        mass_percent,
        dsc_mW_per_mg,
        purge_flow_mL_per_min,
        protective_flow_mL_per_min,
        sensitivity_uV_per_mW,
        segment,
        nominal_heating_rate_K_per_min,
        nominal_hold_temperature_celsius,
        atmosphere,
        initial_mass_mg,
        acquisition_datetime,
        raw_export_file,
        legacy_alias,
        replicate,
        notes,
        import_warnings,
    )
        experiment = new(
            id,
            role,
            kind,
            source_file,
            source_variable,
            source_sha256,
            sample_id,
            composition,
            temperature_K,
            time_min,
            mass_percent,
            dsc_mW_per_mg,
            purge_flow_mL_per_min,
            protective_flow_mL_per_min,
            sensitivity_uV_per_mW,
            segment,
            nominal_heating_rate_K_per_min,
            nominal_hold_temperature_celsius,
            atmosphere,
            initial_mass_mg,
            acquisition_datetime,
            raw_export_file,
            legacy_alias,
            replicate,
            notes,
            import_warnings,
        )
        issues = validate_experiment(experiment)
        errors = filter(issue -> issue.severity == :error, issues)
        isempty(errors) || throw(ExperimentValidationError(experiment.id, errors))
        for issue in filter(issue -> issue.severity == :warning, issues)
            warning = "$(issue.code): $(issue.message)"
            warning in experiment.import_warnings ||
                push!(experiment.import_warnings, warning)
        end
        return experiment
    end
end

function _optional_float_vector(values)
    return isnothing(values) ? nothing : Float64.(collect(values))
end

"""
Construct and validate a raw experiment from named fields.
"""
function Experiment(;
    id::AbstractString,
    role::Symbol,
    kind::Symbol,
    source_file::AbstractString,
    source_variable::AbstractString,
    source_sha256::AbstractString,
    sample_id::AbstractString,
    composition::Composition,
    temperature_K,
    time_min,
    mass_percent,
    dsc_mW_per_mg=nothing,
    purge_flow_mL_per_min=nothing,
    protective_flow_mL_per_min=nothing,
    sensitivity_uV_per_mW=nothing,
    segment=nothing,
    nominal_heating_rate_K_per_min::Real,
    nominal_hold_temperature_celsius=nothing,
    atmosphere::AbstractString,
    initial_mass_mg=missing,
    acquisition_datetime=missing,
    raw_export_file=nothing,
    legacy_alias=nothing,
    replicate=missing,
    notes=String[],
    import_warnings=String[],
)
    hold_temperature = if isnothing(nominal_hold_temperature_celsius)
        nothing
    else
        Float64(nominal_hold_temperature_celsius)
    end
    sample_mass = ismissing(initial_mass_mg) ? missing : Float64(initial_mass_mg)
    acquired = ismissing(acquisition_datetime) ? missing : DateTime(acquisition_datetime)
    replicate_label = ismissing(replicate) ? missing : String(replicate)
    segments = if isnothing(segment)
        nothing
    else
        Union{Missing,Int}[value for value in segment]
    end

    return Experiment(
        Val(:validated),
        String(id),
        role,
        kind,
        String(source_file),
        String(source_variable),
        lowercase(String(source_sha256)),
        String(sample_id),
        composition,
        Float64.(collect(temperature_K)),
        Float64.(collect(time_min)),
        Float64.(collect(mass_percent)),
        _optional_float_vector(dsc_mW_per_mg),
        _optional_float_vector(purge_flow_mL_per_min),
        _optional_float_vector(protective_flow_mL_per_min),
        _optional_float_vector(sensitivity_uV_per_mW),
        segments,
        Float64(nominal_heating_rate_K_per_min),
        hold_temperature,
        String(atmosphere),
        sample_mass,
        acquired,
        isnothing(raw_export_file) ? nothing : String(raw_export_file),
        isnothing(legacy_alias) ? nothing : String(legacy_alias),
        replicate_label,
        String.(collect(notes)),
        String.(collect(import_warnings)),
    )
end

"""
Boolean mask of rows with finite temperature, time, and mass values.
"""
function valid_row_mask(experiment::Experiment)
    return isfinite.(experiment.temperature_K) .& isfinite.(experiment.time_min) .&
           isfinite.(experiment.mass_percent)
end

"""
Number of consecutive invalid core rows at the end of an experiment.
"""
function trailing_invalid_row_count(experiment::Experiment)
    mask = valid_row_mask(experiment)
    count = 0
    for valid in Iterators.reverse(mask)
        valid && break
        count += 1
    end
    return count
end

function _channel_length_issues!(issues, experiment, expected_length)
    channels = (
        dsc_mW_per_mg=experiment.dsc_mW_per_mg,
        purge_flow_mL_per_min=experiment.purge_flow_mL_per_min,
        protective_flow_mL_per_min=experiment.protective_flow_mL_per_min,
        sensitivity_uV_per_mW=experiment.sensitivity_uV_per_mW,
        segment=experiment.segment,
    )
    for (name, values) in pairs(channels)
        if !isnothing(values) && length(values) != expected_length
            push!(
                issues,
                ValidationIssue(
                    :error,
                    :channel_length_mismatch,
                    "$name has $(length(values)) rows; expected $expected_length",
                ),
            )
        end
    end
    return issues
end

"""
    validate_experiment(experiment) -> Vector{ValidationIssue}

Validate identity, composition, array lengths, finite-row availability, ordering, and basic
metadata. Raw non-monotonicity and invalid acquisition rows are reported, not removed.
"""
function validate_experiment(experiment::Experiment)
    issues = ValidationIssue[]
    isempty(strip(experiment.id)) &&
        push!(issues, ValidationIssue(:error, :missing_id, "id cannot be empty"))
    experiment.role in (:calibration, :validation) || push!(
        issues,
        ValidationIssue(:error, :invalid_role, "role must be calibration or validation"),
    )
    experiment.kind in (:dynamic, :ramp_hold) || push!(
        issues,
        ValidationIssue(:error, :invalid_kind, "kind must be dynamic or ramp_hold"),
    )
    isempty(strip(experiment.source_file)) && push!(
        issues,
        ValidationIssue(:error, :missing_source_file, "source_file cannot be empty"),
    )
    isempty(strip(experiment.source_variable)) && push!(
        issues,
        ValidationIssue(
            :error, :missing_source_variable, "source_variable cannot be empty"
        ),
    )
    occursin(r"^[0-9a-f]{64}$", experiment.source_sha256) || push!(
        issues,
        ValidationIssue(
            :error, :invalid_source_hash, "source_sha256 must contain 64 hex digits"
        ),
    )
    isfinite(experiment.nominal_heating_rate_K_per_min) &&
    experiment.nominal_heating_rate_K_per_min > 0 || push!(
        issues,
        ValidationIssue(
            :error,
            :invalid_nominal_heating_rate,
            "nominal heating rate must be finite and positive",
        ),
    )
    if experiment.kind == :ramp_hold &&
        isnothing(experiment.nominal_hold_temperature_celsius)
        push!(
            issues,
            ValidationIssue(
                :error,
                :missing_hold_temperature,
                "ramp_hold experiments require a nominal hold temperature",
            ),
        )
    end

    row_count = length(experiment.temperature_K)
    row_count == length(experiment.time_min) == length(experiment.mass_percent) || push!(
        issues,
        ValidationIssue(
            :error,
            :core_length_mismatch,
            "temperature, time, and mass channels must have equal lengths",
        ),
    )
    row_count > 1 || push!(
        issues,
        ValidationIssue(:error, :insufficient_rows, "at least two rows are required"),
    )
    _channel_length_issues!(issues, experiment, row_count)

    if row_count == length(experiment.time_min) == length(experiment.mass_percent)
        mask = valid_row_mask(experiment)
        valid_indices = findall(mask)
        length(valid_indices) >= 2 || push!(
            issues,
            ValidationIssue(
                :error,
                :insufficient_finite_rows,
                "at least two finite core rows are required",
            ),
        )
        invalid_count = count(!, mask)
        trailing_count = trailing_invalid_row_count(experiment)
        if trailing_count > 0
            push!(
                issues,
                ValidationIssue(
                    :warning,
                    :trailing_invalid_rows,
                    "$trailing_count trailing row(s) have non-finite temperature, time, or mass",
                ),
            )
        end
        internal_count = invalid_count - trailing_count
        if internal_count > 0
            push!(
                issues,
                ValidationIssue(
                    :warning,
                    :internal_invalid_rows,
                    "$internal_count non-trailing row(s) have non-finite temperature, time, or mass",
                ),
            )
        end
        if length(valid_indices) >= 2
            times = experiment.time_min[valid_indices]
            time_nonincreasing = count(delta -> delta <= 0, diff(times))
            time_nonincreasing > 0 && push!(
                issues,
                ValidationIssue(
                    :warning,
                    :nonincreasing_time,
                    "$time_nonincreasing finite time step(s) are zero or negative",
                ),
            )
            temperatures = experiment.temperature_K[valid_indices]
            temperature_nonincreasing = count(delta -> delta <= 0, diff(temperatures))
            temperature_nonincreasing > 0 && push!(
                issues,
                ValidationIssue(
                    :info,
                    :nonincreasing_temperature,
                    "$temperature_nonincreasing finite temperature step(s) are zero or negative",
                ),
            )
        end
    end
    return issues
end
