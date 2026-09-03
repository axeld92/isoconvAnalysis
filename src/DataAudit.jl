"""
Deterministic ingestion diagnostics for one raw experiment.
"""
struct ExperimentAudit
    id::String
    role::Symbol
    kind::Symbol
    source_variable::String
    waste_tire_dry_mass_percent::Float64
    nominal_heating_rate_K_per_min::Float64
    measured_heating_rate_K_per_min::Float64
    nominal_hold_temperature_celsius::Union{Nothing,Float64}
    row_count::Int
    finite_core_row_count::Int
    invalid_core_row_count::Int
    trailing_invalid_row_count::Int
    time_nonincreasing_step_count::Int
    temperature_nonincreasing_step_count::Int
    median_sampling_interval_min::Float64
    median_purge_flow_mL_per_min::Union{Nothing,Float64}
    median_protective_flow_mL_per_min::Union{Nothing,Float64}
    minimum_temperature_celsius::Float64
    maximum_temperature_celsius::Float64
    minimum_time_min::Float64
    maximum_time_min::Float64
    minimum_mass_percent::Float64
    maximum_mass_percent::Float64
    segment_values::Vector{Int}
    warnings::Vector{String}
end

function _finite_median(values)
    isnothing(values) && return nothing
    finite_values = filter(isfinite, values)
    return isempty(finite_values) ? nothing : median(finite_values)
end

function _ols_slope(x::AbstractVector, y::AbstractVector)
    length(x) >= 2 || return NaN
    centered_x = x .- mean(x)
    denominator = sum(abs2, centered_x)
    denominator > 0 || return NaN
    return sum(centered_x .* (y .- mean(y))) / denominator
end

function _heating_rate_indices(experiment::Experiment, valid_mask::BitVector)
    if experiment.kind == :ramp_hold && !isnothing(experiment.segment)
        return findall(
            index -> valid_mask[index] && experiment.segment[index] === 3,
            eachindex(valid_mask),
        )
    end
    temperature_celsius = experiment.temperature_K .- 273.15
    return findall(
        index -> valid_mask[index] && 100.0 <= temperature_celsius[index] <= 700.0,
        eachindex(valid_mask),
    )
end

"""
Compute import, range, ordering, sampling, segment, and measured-rate diagnostics.
"""
function audit_experiment(experiment::Experiment)
    mask = BitVector(valid_row_mask(experiment))
    valid_indices = findall(mask)
    temperatures_celsius = experiment.temperature_K[valid_indices] .- 273.15
    times = experiment.time_min[valid_indices]
    masses = experiment.mass_percent[valid_indices]
    time_steps = diff(times)
    temperature_steps = diff(temperatures_celsius)
    positive_time_steps = filter(>(0), time_steps)
    rate_indices = _heating_rate_indices(experiment, mask)
    measured_rate = _ols_slope(
        experiment.time_min[rate_indices], experiment.temperature_K[rate_indices]
    )
    segments = if isnothing(experiment.segment)
        Int[]
    else
        sort!(unique(Int[value for value in experiment.segment if !ismissing(value)]))
    end

    return ExperimentAudit(
        experiment.id,
        experiment.role,
        experiment.kind,
        experiment.source_variable,
        100 * experiment.composition.waste_tire,
        experiment.nominal_heating_rate_K_per_min,
        measured_rate,
        experiment.nominal_hold_temperature_celsius,
        length(experiment.temperature_K),
        length(valid_indices),
        count(!, mask),
        trailing_invalid_row_count(experiment),
        count(delta -> delta <= 0, time_steps),
        count(delta -> delta <= 0, temperature_steps),
        isempty(positive_time_steps) ? NaN : median(positive_time_steps),
        _finite_median(experiment.purge_flow_mL_per_min),
        _finite_median(experiment.protective_flow_mL_per_min),
        minimum(temperatures_celsius),
        maximum(temperatures_celsius),
        minimum(times),
        maximum(times),
        minimum(masses),
        maximum(masses),
        segments,
        copy(experiment.import_warnings),
    )
end

"""
Load and audit every experiment selected by the catalog.
"""
function audit_catalog(catalog::DatasetCatalog)
    return audit_experiment.(load_experiments(catalog))
end

_report_float(value; digits=4) = isfinite(value) ? @sprintf("%.*f", digits, value) : "n/a"

function _flow_summary(audits, kind)
    pairs = unique([
        (audit.median_purge_flow_mL_per_min, audit.median_protective_flow_mL_per_min) for
        audit in audits if audit.kind == kind
    ])
    labels = map(pairs) do (purge, protective)
        purge_label = isnothing(purge) ? "n/a" : _report_float(purge; digits=1)
        protective_label =
            isnothing(protective) ? "n/a" : _report_float(protective; digits=1)
        return "$purge_label/$protective_label mL/min"
    end
    return join(labels, ", ")
end

function _relative_to_project(catalog::DatasetCatalog, path::String)
    project_root = normpath(joinpath(dirname(catalog.source_path), ".."))
    return relpath(path, project_root)
end

function _source_inventory(catalog::DatasetCatalog)
    return map(catalog.sources) do source
        present = list_mat_variables(source.source_file)
        selected = getfield.(source.experiments, :source_variable)
        return (
            id=source.id,
            file=_relative_to_project(catalog, source.source_file),
            sha256=source.expected_sha256,
            present=present,
            selected=selected,
            missing=sort!(collect(setdiff(Set(selected), Set(present)))),
            unselected=sort!(collect(setdiff(Set(present), Set(selected)))),
        )
    end
end

function _audit_toml(catalog::DatasetCatalog, audits::Vector{ExperimentAudit})
    sources = _source_inventory(catalog)
    experiment_tables = map(audits) do audit
        table = Dict{String,Any}(
            "id" => audit.id,
            "role" => String(audit.role),
            "kind" => String(audit.kind),
            "source_variable" => audit.source_variable,
            "waste_tire_dry_mass_percent" => audit.waste_tire_dry_mass_percent,
            "nominal_heating_rate_K_per_min" => audit.nominal_heating_rate_K_per_min,
            "measured_heating_rate_K_per_min" => audit.measured_heating_rate_K_per_min,
            "row_count" => audit.row_count,
            "finite_core_row_count" => audit.finite_core_row_count,
            "invalid_core_row_count" => audit.invalid_core_row_count,
            "trailing_invalid_row_count" => audit.trailing_invalid_row_count,
            "time_nonincreasing_step_count" => audit.time_nonincreasing_step_count,
            "temperature_nonincreasing_step_count" =>
                audit.temperature_nonincreasing_step_count,
            "median_sampling_interval_min" => audit.median_sampling_interval_min,
            "minimum_temperature_celsius" => audit.minimum_temperature_celsius,
            "maximum_temperature_celsius" => audit.maximum_temperature_celsius,
            "minimum_time_min" => audit.minimum_time_min,
            "maximum_time_min" => audit.maximum_time_min,
            "minimum_mass_percent" => audit.minimum_mass_percent,
            "maximum_mass_percent" => audit.maximum_mass_percent,
            "segment_values" => audit.segment_values,
            "warnings" => audit.warnings,
        )
        if !isnothing(audit.nominal_hold_temperature_celsius)
            table["nominal_hold_temperature_celsius"] =
                audit.nominal_hold_temperature_celsius
        end
        if !isnothing(audit.median_purge_flow_mL_per_min)
            table["median_purge_flow_mL_per_min"] = audit.median_purge_flow_mL_per_min
        end
        if !isnothing(audit.median_protective_flow_mL_per_min)
            table["median_protective_flow_mL_per_min"] =
                audit.median_protective_flow_mL_per_min
        end
        return table
    end
    source_tables = [
        Dict{String,Any}(
            "id" => source.id,
            "file" => source.file,
            "sha256" => source.sha256,
            "present_variables" => source.present,
            "selected_variables" => source.selected,
            "missing_variables" => source.missing,
            "unselected_variables" => source.unselected,
        ) for source in sources
    ]
    return Dict{String,Any}(
        "schema_version" => 1,
        "catalog_file" => _relative_to_project(catalog, catalog.source_path),
        "summary" => Dict{String,Any}(
            "source_count" => length(catalog.sources),
            "experiment_count" => length(audits),
            "calibration_experiment_count" =>
                count(audit -> audit.role == :calibration, audits),
            "validation_experiment_count" =>
                count(audit -> audit.role == :validation, audits),
            "total_row_count" => sum(getfield.(audits, :row_count)),
            "invalid_core_row_count" => sum(getfield.(audits, :invalid_core_row_count)),
            "missing_selected_variable_count" =>
                sum(length(source.missing) for source in sources),
        ),
        "sources" => source_tables,
        "experiments" => experiment_tables,
    )
end

function _write_markdown_audit(
    io::IO, catalog::DatasetCatalog, audits::Vector{ExperimentAudit}
)
    inventories = _source_inventory(catalog)
    println(io, "# M2 raw-data audit")
    println(io)
    println(io, "Generated by `scripts/audit_legacy_data.jl` from `config/datasets.toml`.")
    println(io, "No acquisition row was removed or modified during this audit.")
    println(io)
    println(io, "## Summary")
    println(io)
    println(io, "- Sources: $(length(inventories)) MAT v5 files")
    println(io, "- Selected experiments: $(length(audits))")
    println(
        io,
        "- Dynamic calibration experiments: $(count(a -> a.role == :calibration, audits))",
    )
    println(
        io,
        "- Ramp-and-hold validation experiments: $(count(a -> a.role == :validation, audits))",
    )
    println(io, "- Total stored rows: $(sum(a.row_count for a in audits))")
    println(io, "- Invalid core rows: $(sum(a.invalid_core_row_count for a in audits))")
    println(
        io,
        "- Missing selected MAT variables: $(sum(length(s.missing) for s in inventories))",
    )
    println(io)
    println(
        io,
        "All $(length(audits)) selected experiments loaded and passed structural validation.",
    )
    println(
        io,
        "The raw arrays retain $(sum(a.trailing_invalid_row_count for a in audits)) trailing invalid core row(s); none were removed.",
    )
    println(
        io,
        "Finite time channels contain $(sum(a.time_nonincreasing_step_count for a in audits)) zero or negative step(s).",
    )
    println(io)
    println(io, "## Source inventory")
    println(io)
    println(io, "| Source | File | Present | Selected | Missing | Unselected | SHA-256 |")
    println(io, "|---|---|---:|---:|---:|---|---|")
    for source in inventories
        unselected = isempty(source.unselected) ? "—" : join(source.unselected, ", ")
        println(
            io,
            "| `$(source.id)` | `$(source.file)` | $(length(source.present)) | $(length(source.selected)) | $(length(source.missing)) | $unselected | `$(source.sha256)` |",
        )
    end
    println(io)
    if any("ExpDatdrifttest2" in source.unselected for source in inventories)
        println(
            io,
            "`ExpDatdrifttest2` is present in the dynamic MAT file but deliberately unselected:",
        )
        println(
            io,
            "it is an instrument drift test, not one of the 20 composition/rate conditions.",
        )
        println(io)
    end
    println(io, "## Dynamic calibration runs")
    println(io)
    println(
        io,
        "| ID | Variable | Tire wt% | Nominal β | Measured β | Rows | Invalid tail | T range (°C) | Mass range (%) | Δt median (s) | Segments |",
    )
    println(io, "|---|---|---:|---:|---:|---:|---:|---|---|---:|---|")
    for audit in filter(audit -> audit.kind == :dynamic, audits)
        segments =
            isempty(audit.segment_values) ? "absent" : join(audit.segment_values, ",")
        println(
            io,
            "| `$(audit.id)` | `$(audit.source_variable)` | $(_report_float(audit.waste_tire_dry_mass_percent; digits=0)) | $(_report_float(audit.nominal_heating_rate_K_per_min; digits=1)) | $(_report_float(audit.measured_heating_rate_K_per_min)) | $(audit.row_count) | $(audit.trailing_invalid_row_count) | $(_report_float(audit.minimum_temperature_celsius; digits=2))–$(_report_float(audit.maximum_temperature_celsius; digits=2)) | $(_report_float(audit.minimum_mass_percent; digits=2))–$(_report_float(audit.maximum_mass_percent; digits=2)) | $(_report_float(60 * audit.median_sampling_interval_min; digits=3)) | $segments |",
        )
    end
    println(io)
    println(
        io,
        "Measured dynamic rates are OLS slopes of temperature versus time over 100–700 °C;",
    )
    println(io, "they are audit diagnostics, not preprocessing outputs.")
    println(io)
    println(io, "## Ramp-and-hold validation runs")
    println(io)
    println(
        io,
        "| ID | Variable | Tire wt% | Hold (°C) | Ramp β | Rows | Invalid tail | T range (°C) | Mass range (%) | Δt median (s) | Segments |",
    )
    println(io, "|---|---|---:|---:|---:|---:|---:|---|---|---:|---|")
    for audit in filter(audit -> audit.kind == :ramp_hold, audits)
        segments = join(audit.segment_values, ",")
        println(
            io,
            "| `$(audit.id)` | `$(audit.source_variable)` | $(_report_float(audit.waste_tire_dry_mass_percent; digits=0)) | $(_report_float(something(audit.nominal_hold_temperature_celsius); digits=0)) | $(_report_float(audit.measured_heating_rate_K_per_min)) | $(audit.row_count) | $(audit.trailing_invalid_row_count) | $(_report_float(audit.minimum_temperature_celsius; digits=2))–$(_report_float(audit.maximum_temperature_celsius; digits=2)) | $(_report_float(audit.minimum_mass_percent; digits=2))–$(_report_float(audit.maximum_mass_percent; digits=2)) | $(_report_float(60 * audit.median_sampling_interval_min; digits=3)) | $segments |",
        )
    end
    println(io)
    println(
        io,
        "Ramp rates are OLS slopes within recorded segment 3. Segment 4 is the hold and is",
    )
    println(io, "not included in the rate estimate.")
    println(io)
    println(io, "## Findings carried into M3")
    println(io)
    println(io, "- The importer preserves all trailing invalid rows in selected matrices.")
    low_rate_endocarp_ids = Set((
        "dynamic_wt000_rate05", "dynamic_wt000_rate10", "dynamic_wt000_rate15"
    ))
    if issubset(low_rate_endocarp_ids, Set(getfield.(audits, :id)))
        println(
            io,
            "- The 0 wt% 5, 10, and 15 K/min matrices omit the optional segment channel.",
        )
        println(
            io,
            "  They begin near time 30 min and mass 94%, so import must not infer an initial",
        )
        println(io, "  100% mass point or re-zero time.")
    end
    any(audit -> audit.kind == :dynamic, audits) && println(
        io,
        "- Recorded median purge/protective flows are $(_flow_summary(audits, :dynamic)) for dynamic runs.",
    )
    any(audit -> audit.kind == :ramp_hold, audits) && println(
        io,
        "- Recorded median purge/protective flows are $(_flow_summary(audits, :ramp_hold)) for ramp-and-hold runs.",
    )
    if sum(audit.temperature_nonincreasing_step_count for audit in audits) > 0
        println(
            io,
            "- Temperature is not strictly increasing because of conditioning/hold segments and",
        )
        println(
            io,
            "  instrument-scale fluctuations; the machine-readable audit records exact counts.",
        )
    end
    if any(source.id == "ramp_hold_2021" for source in catalog.sources)
        println(
            io,
            "- Ramp-and-hold auxiliary channel labels remain provisional as documented in the data",
        )
        println(
            io,
            "  dictionary. Temperature, time, mass, and segment columns are sufficient for M2.",
        )
        println(
            io,
            "- No replicate uncertainty can be estimated from the retained dynamic design.",
        )
    end
end

"""
    write_audit_reports(catalog, audits; markdown_path, toml_path)

Write deterministic human-readable Markdown and machine-readable TOML audit evidence.
"""
function write_audit_reports(
    catalog::DatasetCatalog,
    audits::Vector{ExperimentAudit};
    markdown_path::AbstractString,
    toml_path::AbstractString,
)
    markdown_output = abspath(markdown_path)
    toml_output = abspath(toml_path)
    mkpath(dirname(markdown_output))
    mkpath(dirname(toml_output))
    open(markdown_output, "w") do io
        return _write_markdown_audit(io, catalog, audits)
    end
    open(toml_output, "w") do io
        return TOML.print(io, _audit_toml(catalog, audits); sorted=true)
    end
    return (markdown=markdown_output, toml=toml_output)
end
