const _SUPPORTED_BOUNDARY_COLUMNS = Set((
    :temperature_C,
    :time_min,
    :dsc_mW_per_mg,
    :mass_percent,
    :purge_flow_mL_per_min,
    :protective_flow_mL_per_min,
    :sensitivity_uV_per_mW,
    :segment,
))

function _dataset_error(path::AbstractString, message::AbstractString)
    return throw(DatasetConfigError(String(path), String(message)))
end

function _dataset_required(table::AbstractDict, key::String, context::String, path::String)
    haskey(table, key) || _dataset_error(path, "missing $context.$key")
    return table[key]
end

function _dataset_string(table::AbstractDict, key::String, context::String, path::String)
    value = _dataset_required(table, key, context, path)
    value isa AbstractString || _dataset_error(path, "$context.$key must be a string")
    isempty(strip(value)) && _dataset_error(path, "$context.$key cannot be empty")
    return String(value)
end

function _dataset_number(table::AbstractDict, key::String, context::String, path::String)
    value = _dataset_required(table, key, context, path)
    value isa Real && !(value isa Bool) ||
        _dataset_error(path, "$context.$key must be a number")
    isfinite(value) || _dataset_error(path, "$context.$key must be finite")
    return Float64(value)
end

function _dataset_strings(
    table::AbstractDict, key::String, context::String, path::String; required::Bool=true
)
    if !haskey(table, key)
        required && _dataset_error(path, "missing $context.$key")
        return String[]
    end
    values = table[key]
    values isa AbstractVector || _dataset_error(path, "$context.$key must be an array")
    all(value -> value isa AbstractString, values) ||
        _dataset_error(path, "$context.$key must contain only strings")
    return String.(values)
end

function _dataset_optional_string(
    table::AbstractDict, key::String, context::String, path::String
)
    haskey(table, key) || return nothing
    return _dataset_string(table, key, context, path)
end

function _dataset_optional_number(
    table::AbstractDict, key::String, context::String, path::String
)
    haskey(table, key) || return nothing
    return _dataset_number(table, key, context, path)
end

function _resolve_catalog_path(catalog_path::String, value::String)
    return if isabspath(value)
        normpath(value)
    else
        normpath(joinpath(dirname(catalog_path), value))
    end
end

function _parse_experiment_spec(
    table::AbstractDict, source_context::String, catalog_path::String
)
    context = "$source_context.experiments"
    id = _dataset_string(table, "id", context, catalog_path)
    variable = _dataset_string(table, "variable", context, catalog_path)
    sample_id = _dataset_string(table, "sample_id", context, catalog_path)
    waste_tire_percent = _dataset_number(
        table, "waste_tire_dry_mass_percent", context, catalog_path
    )
    0.0 <= waste_tire_percent <= 100.0 || _dataset_error(
        catalog_path, "$context.waste_tire_dry_mass_percent must lie in [0, 100]"
    )
    nominal_rate = _dataset_number(
        table, "nominal_heating_rate_K_per_min", context, catalog_path
    )
    nominal_rate > 0 || _dataset_error(
        catalog_path, "$context.nominal_heating_rate_K_per_min must be positive"
    )

    hold_temperature = _dataset_optional_number(
        table, "nominal_hold_temperature_celsius", context, catalog_path
    )
    initial_mass = _dataset_optional_number(table, "initial_mass_mg", context, catalog_path)
    acquired_text = _dataset_optional_string(
        table, "acquisition_datetime", context, catalog_path
    )
    acquired_at = if isnothing(acquired_text)
        missing
    else
        try
            DateTime(acquired_text)
        catch
            _dataset_error(
                catalog_path,
                "$context.acquisition_datetime must use ISO format yyyy-mm-ddTHH:MM:SS",
            )
        end
    end
    raw_export = _dataset_optional_string(table, "raw_export_file", context, catalog_path)
    raw_export_path =
        isnothing(raw_export) ? nothing : _resolve_catalog_path(catalog_path, raw_export)
    if !isnothing(raw_export_path) && !isfile(raw_export_path)
        _dataset_error(
            catalog_path, "$context.raw_export_file does not exist: $raw_export_path"
        )
    end
    legacy_alias = _dataset_optional_string(table, "legacy_alias", context, catalog_path)
    replicate = _dataset_optional_string(table, "replicate", context, catalog_path)
    notes = _dataset_strings(table, "notes", context, catalog_path; required=false)

    return ExperimentSpec(
        id,
        variable,
        sample_id,
        Composition(waste_tire_percent / 100),
        nominal_rate,
        hold_temperature,
        isnothing(initial_mass) ? missing : initial_mass,
        acquired_at,
        raw_export_path,
        legacy_alias,
        isnothing(replicate) ? missing : replicate,
        notes,
    )
end

function _parse_source_spec(table::AbstractDict, index::Int, catalog_path::String)
    context = "sources[$index]"
    id = _dataset_string(table, "id", context, catalog_path)
    role = Symbol(_dataset_string(table, "role", context, catalog_path))
    role in (:calibration, :validation) ||
        _dataset_error(catalog_path, "$context.role must be calibration or validation")
    kind = Symbol(_dataset_string(table, "kind", context, catalog_path))
    kind in (:dynamic, :ramp_hold) ||
        _dataset_error(catalog_path, "$context.kind must be dynamic or ramp_hold")
    format = Symbol(_dataset_string(table, "format", context, catalog_path))
    format == Symbol("mat-v5") ||
        _dataset_error(catalog_path, "$context.format must be mat-v5")
    source_file = _resolve_catalog_path(
        catalog_path, _dataset_string(table, "source_file", context, catalog_path)
    )
    isfile(source_file) ||
        _dataset_error(catalog_path, "$context.source_file does not exist: $source_file")
    expected_sha256 = lowercase(_dataset_string(table, "sha256", context, catalog_path))
    occursin(r"^[0-9a-f]{64}$", expected_sha256) ||
        _dataset_error(catalog_path, "$context.sha256 must contain 64 hex digits")
    atmosphere = _dataset_string(table, "atmosphere", context, catalog_path)
    required_columns = Symbol.(
        _dataset_strings(table, "required_columns", context, catalog_path)
    )
    optional_columns = Symbol.(
        _dataset_strings(table, "optional_columns", context, catalog_path; required=false)
    )
    isempty(required_columns) &&
        _dataset_error(catalog_path, "$context.required_columns cannot be empty")
    length(unique(vcat(required_columns, optional_columns))) ==
    length(required_columns) + length(optional_columns) ||
        _dataset_error(catalog_path, "$context column names must be unique")
    all(
        column -> column in _SUPPORTED_BOUNDARY_COLUMNS,
        vcat(required_columns, optional_columns),
    ) || _dataset_error(catalog_path, "$context contains an unsupported boundary column")
    all(column -> column in required_columns, (:temperature_C, :time_min, :mass_percent)) ||
        _dataset_error(
            catalog_path,
            "$context.required_columns must include temperature_C, time_min, and mass_percent",
        )
    notes = _dataset_strings(table, "notes", context, catalog_path; required=false)

    experiment_tables = _dataset_required(table, "experiments", context, catalog_path)
    experiment_tables isa AbstractVector ||
        _dataset_error(catalog_path, "$context.experiments must be an array of tables")
    isempty(experiment_tables) &&
        _dataset_error(catalog_path, "$context.experiments cannot be empty")
    experiments = ExperimentSpec[
        _parse_experiment_spec(experiment, context, catalog_path) for
        experiment in experiment_tables
    ]
    ids = getfield.(experiments, :id)
    length(unique(ids)) == length(ids) ||
        _dataset_error(catalog_path, "$context contains duplicate experiment ids")
    variables = getfield.(experiments, :source_variable)
    length(unique(variables)) == length(variables) ||
        _dataset_error(catalog_path, "$context contains duplicate source variables")
    if kind == :ramp_hold &&
        any(isnothing(spec.nominal_hold_temperature_celsius) for spec in experiments)
        _dataset_error(
            catalog_path,
            "$context ramp_hold experiments must declare nominal_hold_temperature_celsius",
        )
    end

    return DatasetSourceSpec(
        id,
        role,
        kind,
        format,
        source_file,
        expected_sha256,
        atmosphere,
        required_columns,
        optional_columns,
        notes,
        experiments,
    )
end

"""
    load_dataset_catalog(path) -> DatasetCatalog

Parse, validate, and resolve every path in a versioned dataset catalog. The catalog names
each MAT variable explicitly; no workspace or filename-pattern import is performed.
"""
function load_dataset_catalog(path::AbstractString)
    catalog_path = abspath(path)
    isfile(catalog_path) || _dataset_error(catalog_path, "file does not exist")
    data = try
        TOML.parsefile(catalog_path)
    catch error
        error isa TOML.ParserError || rethrow()
        _dataset_error(catalog_path, sprint(showerror, error))
    end
    schema_version = _dataset_required(data, "schema_version", "catalog", catalog_path)
    schema_version isa Integer ||
        _dataset_error(catalog_path, "catalog.schema_version must be an integer")
    schema_version == 1 ||
        _dataset_error(catalog_path, "unsupported schema_version $schema_version")
    source_tables = _dataset_required(data, "sources", "catalog", catalog_path)
    source_tables isa AbstractVector ||
        _dataset_error(catalog_path, "catalog.sources must be an array of tables")
    sources = DatasetSourceSpec[
        _parse_source_spec(source, index, catalog_path) for
        (index, source) in enumerate(source_tables)
    ]
    source_ids = getfield.(sources, :id)
    length(unique(source_ids)) == length(source_ids) ||
        _dataset_error(catalog_path, "source ids must be unique")
    experiment_ids = [
        experiment.id for source in sources for experiment in source.experiments
    ]
    length(unique(experiment_ids)) == length(experiment_ids) ||
        _dataset_error(catalog_path, "experiment ids must be unique across the catalog")
    return DatasetCatalog(catalog_path, Int(schema_version), sources)
end

"""
Calculate the lowercase SHA-256 fingerprint of a source file.
"""
function source_sha256(path::AbstractString)
    isfile(path) || throw(DataImportError(String(path), "", "file does not exist"))
    return open(path, "r") do io
        return bytes2hex(SHA.sha256(io))
    end
end

function _verified_source_hash(source::DatasetSourceSpec)
    actual = source_sha256(source.source_file)
    actual == source.expected_sha256 || throw(
        DataImportError(
            source.source_file,
            "",
            "SHA-256 mismatch: expected $(source.expected_sha256), received $actual",
        ),
    )
    return actual
end

"""
List the variables stored in a MAT file without loading their values.
"""
function list_mat_variables(path::AbstractString)
    isfile(path) || throw(DataImportError(String(path), "", "file does not exist"))
    file = MAT.matopen(path)
    try
        return sort!(String.(collect(keys(file))))
    finally
        close(file)
    end
end

function _read_mat_matrix(file, source_file::String, variable::String)
    variable in keys(file) ||
        throw(DataImportError(source_file, variable, "variable was not found"))
    value = try
        read(file, variable)
    catch error
        throw(DataImportError(source_file, variable, sprint(showerror, error)))
    end
    value isa AbstractMatrix{<:Real} || throw(
        DataImportError(
            source_file,
            variable,
            "expected a real numeric matrix, received $(typeof(value))",
        ),
    )
    return Matrix{Float64}(value)
end

"""
    load_mat_variable(path, variable) -> Matrix{Float64}

Load exactly one named real matrix from a MAT file. Missing or nonmatrix variables fail
explicitly rather than causing the complete MATLAB workspace to be imported.
"""
function load_mat_variable(path::AbstractString, variable::AbstractString)
    source_file = abspath(path)
    isfile(source_file) ||
        throw(DataImportError(source_file, String(variable), "file does not exist"))
    file = MAT.matopen(source_file)
    try
        return _read_mat_matrix(file, source_file, String(variable))
    finally
        close(file)
    end
end

function _matrix_columns(
    source::DatasetSourceSpec, matrix::AbstractMatrix, variable::String
)
    required_count = length(source.required_columns)
    optional_count = length(source.optional_columns)
    column_count = size(matrix, 2)
    required_count <= column_count <= required_count + optional_count || throw(
        DataImportError(
            source.source_file,
            variable,
            "matrix has $column_count columns; expected $required_count to $(required_count + optional_count)",
        ),
    )
    used_optional_count = column_count - required_count
    return vcat(source.required_columns, source.optional_columns[1:used_optional_count])
end

function _matrix_column(matrix, columns, name::Symbol)
    index = findfirst(==(name), columns)
    return isnothing(index) ? nothing : vec(matrix[:, index])
end

function _segments_from_raw(raw_segment, warnings)
    isnothing(raw_segment) && return nothing
    result = Union{Missing,Int}[]
    invalid_count = 0
    for value in raw_segment
        if isfinite(value) && isinteger(value)
            push!(result, Int(value))
        else
            push!(result, missing)
            !isnan(value) && (invalid_count += 1)
        end
    end
    invalid_count > 0 && push!(
        warnings,
        "invalid_segment_values: $invalid_count finite segment value(s) were not integers",
    )
    return result
end

function _experiment_from_matrix(
    source::DatasetSourceSpec,
    spec::ExperimentSpec,
    matrix::Matrix{Float64},
    actual_hash::String,
)
    columns = _matrix_columns(source, matrix, spec.source_variable)
    warnings = String[]
    if :segment in source.optional_columns && !(:segment in columns)
        push!(warnings, "missing_optional_segment: source matrix omits the segment channel")
    end
    temperature_celsius = _matrix_column(matrix, columns, :temperature_C)
    time_min = _matrix_column(matrix, columns, :time_min)
    mass_percent = _matrix_column(matrix, columns, :mass_percent)
    segment = _segments_from_raw(_matrix_column(matrix, columns, :segment), warnings)

    return Experiment(;
        id=spec.id,
        role=source.role,
        kind=source.kind,
        source_file=source.source_file,
        source_variable=spec.source_variable,
        source_sha256=actual_hash,
        sample_id=spec.sample_id,
        composition=spec.composition,
        temperature_K=temperature_celsius .+ 273.15,
        time_min=time_min,
        mass_percent=mass_percent,
        dsc_mW_per_mg=_matrix_column(matrix, columns, :dsc_mW_per_mg),
        purge_flow_mL_per_min=_matrix_column(matrix, columns, :purge_flow_mL_per_min),
        protective_flow_mL_per_min=_matrix_column(
            matrix, columns, :protective_flow_mL_per_min
        ),
        sensitivity_uV_per_mW=_matrix_column(matrix, columns, :sensitivity_uV_per_mW),
        segment=segment,
        nominal_heating_rate_K_per_min=spec.nominal_heating_rate_K_per_min,
        nominal_hold_temperature_celsius=spec.nominal_hold_temperature_celsius,
        atmosphere=source.atmosphere,
        initial_mass_mg=spec.initial_mass_mg,
        acquisition_datetime=spec.acquisition_datetime,
        raw_export_file=spec.raw_export_file,
        legacy_alias=spec.legacy_alias,
        replicate=spec.replicate,
        notes=vcat(source.notes, spec.notes),
        import_warnings=warnings,
    )
end

function _load_source_experiments(source::DatasetSourceSpec, selected_ids)
    selected_specs = if isnothing(selected_ids)
        source.experiments
    else
        filter(spec -> spec.id in selected_ids, source.experiments)
    end
    isempty(selected_specs) && return Experiment[]
    actual_hash = _verified_source_hash(source)
    file = MAT.matopen(source.source_file)
    try
        return Experiment[
            _experiment_from_matrix(
                source,
                spec,
                _read_mat_matrix(file, source.source_file, spec.source_variable),
                actual_hash,
            ) for spec in selected_specs
        ]
    finally
        close(file)
    end
end

"""
Load one explicitly configured experiment by stable catalog id.
"""
function load_experiment(catalog::DatasetCatalog, id::AbstractString)
    selected_id = String(id)
    matches = [
        (source, spec) for source in catalog.sources for
        spec in source.experiments if spec.id == selected_id
    ]
    length(matches) == 1 ||
        throw(DataImportError(catalog.source_path, selected_id, "catalog id was not found"))
    source, _ = only(matches)
    return only(_load_source_experiments(source, Set((selected_id,))))
end

"""
    load_experiments(catalog; roles=nothing, kinds=nothing, ids=nothing)

Load selected experiments in catalog order, verifying each source checksum once. Optional
filters accept collections of symbols (`roles`, `kinds`) or stable string ids (`ids`).
"""
function load_experiments(
    catalog::DatasetCatalog; roles=nothing, kinds=nothing, ids=nothing
)
    role_filter = isnothing(roles) ? nothing : Set(Symbol.(roles))
    kind_filter = isnothing(kinds) ? nothing : Set(Symbol.(kinds))
    id_filter = isnothing(ids) ? nothing : Set(String.(ids))
    loaded = Experiment[]
    for source in catalog.sources
        !isnothing(role_filter) && !(source.role in role_filter) && continue
        !isnothing(kind_filter) && !(source.kind in kind_filter) && continue
        append!(loaded, _load_source_experiments(source, id_filter))
    end
    if !isnothing(id_filter)
        loaded_ids = Set(getfield.(loaded, :id))
        missing_ids = sort!(collect(setdiff(id_filter, loaded_ids)))
        isempty(missing_ids) || throw(
            DataImportError(
                catalog.source_path,
                join(missing_ids, ", "),
                "one or more requested catalog ids were not found under the selected filters",
            ),
        )
    end
    return loaded
end
