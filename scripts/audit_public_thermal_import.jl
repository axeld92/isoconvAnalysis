using IsoconversionalAnalysis
using SHA
using TOML

const EXPECTED_GROUPS = Dict(
    "synthetic_tg" => 13,
    "experimental_tg" => 5,
    "experimental_dsc_dynamic" => 5,
    "experimental_dsc_isothermal" => 7,
)

function require_check(condition, id, description)
    condition || error("$id: $description")
    return nothing
end

function public_input_path(root, relative)
    components = split(relative, '/')
    valid =
        length(components) > 2 &&
        components[1:2] == ["data", "external"] &&
        all(part -> !isempty(part) && part ∉ (".", ".."), components) &&
        !occursin(r"[\\:]", relative)
    valid || error("Invalid local public input path: $relative")
    source = joinpath(root, components...)
    isfile(source) || error(
        "Missing $relative. Restore the pinned files with scripts/fetch_public_benchmarks.py.",
    )
    resolved = realpath(source)
    data_root = realpath(joinpath(root, "data", "external"))
    relative_to_data = relpath(resolved, data_root)
    (isabspath(relative_to_data) || first(splitpath(relative_to_data)) == "..") &&
        error("Public input resolves outside data/external: $relative")
    return resolved
end

function thermal_mapping(manifest, profile, metadata)
    T = ThermalColumnSpec(
        profile["temperature_column"], Symbol(manifest["temperature_unit"])
    )
    t = ThermalColumnSpec(profile["time_column"], Symbol(manifest["time_unit"]))
    sign = haskey(metadata, "EXO") ? parse(Int, metadata["EXO"]) : nothing
    signals = map(profile["signals"]) do signal
        quantity = Symbol(signal["quantity"])
        return ThermalSignalSpec(
            Symbol(signal["id"]),
            ThermalColumnSpec(signal["column"], Symbol(signal["unit"])),
            quantity;
            exotherm_sign=quantity === :specific_heat_flow ? sign : nothing,
        )
    end
    return ThermalColumnMapping(
        T, signals; time=t, segment_column=get(profile, "segment_column", nothing)
    )
end

function numeric_column_digest(values)
    words = htol.(reinterpret(UInt64, values))
    return bytes2hex(sha256(reinterpret(UInt8, words)))
end

function check_signal(run, signal, expected)
    spec = signal.spec
    require_check(
        String(spec.id) == expected["id"] &&
            spec.column.name == expected["column"] &&
            String(spec.column.unit) == expected["unit"] &&
            String(spec.quantity) == expected["quantity"],
        run.id,
        "signal mapping changed for $(spec.id)",
    )
    values = signal.values
    require_check(
        length(values) == length(run.temperature) && all(isfinite, values),
        run.id,
        "incomplete or nonfinite signal $(spec.id)",
    )
    require_check(signal_values(run, spec.id) == values, run.id, "source signal changed")
    evidence = copy(expected)
    evidence["source_range"] = collect(extrema(values))
    if spec.quantity === :mass
        grams = signal_values(run, spec.id; unit=:g)
        require_check(
            all(isapprox.(grams, values ./ 1_000; rtol=1e-12)),
            run.id,
            "mass did not scale from mg to g",
        )
    elseif spec.quantity === :specific_heat_flow
        watts_per_kg = signal_values(run, spec.id; unit=:W_per_kg)
        require_check(
            watts_per_kg == values .* 1_000,
            run.id,
            "specific heat flow scaling or sign changed",
        )
        require_check(
            spec.exotherm_sign == parse(Int, run.metadata["EXO"]),
            run.id,
            "exotherm convention was lost",
        )
        evidence["exotherm_sign"] = spec.exotherm_sign
    end
    return evidence
end

function check_segments(run, expected_segments)
    id = run.id
    require_check(
        length(run.segments) == length(expected_segments), id, "segment count changed"
    )
    diagnostics = thermal_segment_diagnostics(run)
    first_row = 1
    return map(
        enumerate(zip(run.segments, diagnostics, expected_segments))
    ) do (i, (block, diagnostic, expected))
        count_rows = expected["rows"]
        rows = first_row:(first_row + count_rows - 1)
        require_check(
            block.index == i && block.rows == rows && diagnostic.row_count == count_rows,
            id,
            "segment $i row span changed",
        )
        first_row += count_rows
        require_check(
            block.label == get(expected, "label", nothing), id, "segment $i label changed"
        )
        T = view(run.temperature, rows)
        t = view(run.time, rows)
        primary = view(first(run.signals).values, rows)
        require_check(
            all(>(0), diff(t)), id, "time is not strictly increasing in segment $i"
        )
        for (values, key) in
            ((T, "temperature_range"), (t, "time_range"), (primary, "signal_range"))
            require_check(
                collect(extrema(values)) == expected[key], id, "segment $i $key changed"
            )
        end
        decreasing = count(<(0), diff(T))
        require_check(
            decreasing == expected["temperature_decrease_steps"],
            id,
            "segment $i temperature reversal count changed",
        )
        slope = diagnostic.slope_K_per_min
        require_check(
            slope !== nothing && isapprox(slope, expected["slope"]; rtol=1e-10, atol=1e-11),
            id,
            "segment $i observed temperature/time slope changed",
        )
        evidence = Dict{String,Any}(
            "index" => i,
            "rows" => count_rows,
            "time_range_min" => collect(extrema(t)),
            "temperature_range_C" => collect(extrema(T)),
            "primary_signal_range" => collect(extrema(primary)),
            "temperature_decrease_steps" => decreasing,
            "duration_min" => diagnostic.duration_min,
            "slope_K_per_min" => slope,
        )
        block.label === nothing || (evidence["label"] = block.label)
        return evidence
    end
end

function check_selection(run, expected, segments)
    id = run.id
    index = expected["selected_segment"]
    block = run.segments[index]
    selection = select_thermal_segment(run, index)
    rows = block.rows
    selected_block = only(selection.segments)
    require_check(
        selected_block.index == index &&
            selected_block.label == block.label &&
            selected_block.rows == 1:length(rows),
        id,
        "selected segment identity changed",
    )
    require_check(
        selection.temperature == run.temperature[rows] &&
            selection.time == run.time[rows] &&
            selection.metadata == run.metadata &&
            selection.provenance.sha256 == run.provenance.sha256 &&
            selection.provenance.source_lines == run.provenance.source_lines[rows],
        id,
        "selected segment changed observations or provenance",
    )
    require_check(
        length(selection.signals) == length(run.signals) && all(
            s.spec == r.spec && s.values == r.values[rows] for
            (s, r) in zip(selection.signals, run.signals)
        ),
        id,
        "selected segment changed a recorded signal",
    )

    observed = segments[index]
    result = Dict{String,Any}(
        "index" => index,
        "role" => expected["selected_role"],
        "source_line_range" => collect(extrema(selection.provenance.source_lines)),
    )
    if haskey(expected, "nominal_rate_K_per_min")
        nominal = expected["nominal_rate_K_per_min"]
        require_check(
            isapprox(observed["slope_K_per_min"], nominal; rtol=0.01),
            id,
            "selected first-heating slope does not match the instrument program",
        )
        result["nominal_rate_K_per_min"] = nominal
    elseif haskey(expected, "nominal_temperature_C")
        nominal = expected["nominal_temperature_C"]
        require_check(
            all(T -> abs(T - nominal) <= 2, observed["temperature_range_C"]) &&
                abs(observed["slope_K_per_min"]) < 1e-3,
            id,
            "selected hold does not match the instrument program",
        )
        result["nominal_temperature_C"] = nominal
        program = run.metadata["SEG. $(block.label)"]
        hours, minutes = parse.(Int, split(split(program, '/')[2], ':'))
        planned_minutes = 60 * hours + minutes
        result["planned_duration_min"] = planned_minutes
        result["recorded_duration_shortfall_min"] =
            planned_minutes - observed["duration_min"]
    else
        require_check(
            observed["slope_K_per_min"] > 0, id, "selected TG ramp has no positive slope"
        )
    end
    return result
end

function audit_run(root, manifest, expected)
    id = expected["id"]
    profile = manifest["profiles"][expected["profile"]]
    metadata = get(expected, "metadata", Dict{String,String}())
    mapping = thermal_mapping(manifest, profile, metadata)
    path = public_input_path(root, expected["path"])
    require_check(filesize(path) == expected["bytes"], id, "raw byte count changed")
    run = read_thermal_delimited(
        path,
        mapping;
        id=id,
        material_id=expected["material"],
        delimiter=only(profile["delimiter"]),
        encoding=Symbol(profile["encoding"]),
        header_prefix=get(profile, "header_prefix", nothing),
        expected_sha256=expected["sha256"],
    )
    require_check(
        run.id == id && run.material_id == expected["material"], id, "run identity changed"
    )
    require_check(
        length(run.temperature) == expected["rows"] && length(run.time) == expected["rows"],
        id,
        "row count changed",
    )
    require_check(
        all(isfinite, run.temperature) && all(isfinite, run.time), id, "nonfinite axes"
    )
    require_check(
        run.temperature_unit === :C && temperature_K(run) == run.temperature .+ 273.15,
        id,
        "Celsius to Kelvin conversion changed",
    )
    require_check(
        run.time_unit === :min && time_min(run) == run.time, id, "recorded clock changed"
    )
    provenance = run.provenance
    lines = provenance.source_lines
    expected_first, expected_last = expected["source_data_line_range"]
    require_check(
        provenance.source_path == path &&
            provenance.sha256 == expected["sha256"] &&
            lines == (expected_first:expected_last) &&
            1 < expected_first <= expected_last <= expected["source_lines"],
        id,
        "source path, hash or physical line trace changed",
    )
    require_check(
        isempty(provenance.unmapped_columns), id, "unexpected unmapped source columns"
    )
    require_check(
        length(run.metadata) == expected["source_metadata_count"],
        id,
        "metadata count changed",
    )
    for (key, value) in metadata
        require_check(get(run.metadata, key, nothing) == value, id, "metadata $key changed")
    end
    require_check(
        length(run.signals) == length(profile["signals"]), id, "signal count changed"
    )
    signals = [
        check_signal(run, signal, spec) for
        (signal, spec) in zip(run.signals, profile["signals"])
    ]
    column_digests = Dict(
        "temperature" => numeric_column_digest(run.temperature),
        "time" => numeric_column_digest(run.time),
    )
    for signal in run.signals
        column_digests[String(signal.spec.id)] = numeric_column_digest(signal.values)
    end
    for (column, digest) in expected["numeric_column_sha256"]
        require_check(
            get(column_digests, column, nothing) == digest,
            id,
            "numeric column $column differs from the independent reference",
        )
    end
    require_check(
        length(column_digests) == length(expected["numeric_column_sha256"]),
        id,
        "missing independent numeric column digest",
    )
    segments = check_segments(run, expected["segments"])
    selection = check_selection(run, expected, segments)
    requirements = if expected["selected_role"] == "isothermal_hold"
        [
            "baseline",
            "total_reaction_enthalpy",
            "residual_reaction_enthalpy",
            "initial_conversion",
        ]
    elseif expected["selected_role"] == "first_heating"
        ["baseline", "total_reaction_enthalpy"]
    else
        ["reaction_interval", "initial_mass", "terminal_mass"]
    end

    return Dict{String,Any}(
        "id" => id,
        "material_id" => run.material_id,
        "group" => expected["group"],
        "source_path" => expected["path"],
        "source_sha256" => provenance.sha256,
        "bytes" => filesize(path),
        "rows" => expected["rows"],
        "source_data_line_range" => collect(extrema(lines)),
        "source_metadata_count" => length(run.metadata),
        "numeric_column_sha256" => column_digests,
        "signals" => signals,
        "segments" => segments,
        "selection" => selection,
        "verified_metadata" => metadata,
        "raw_import_verified" => true,
        "kinetic_analysis_ready" => false,
        "conversion_requirements" => requirements,
    )
end

function main()
    root = normpath(joinpath(@__DIR__, ".."))
    manifest_path = "config/public_thermal_import.toml"
    manifest = TOML.parsefile(joinpath(root, manifest_path))
    manifest["schema_version"] == 1 || error("Unsupported public import manifest version")
    runs = manifest["runs"]
    allunique(run["id"] for run in runs) || error("Duplicate public run identifiers")
    allunique(run["path"] for run in runs) || error("Duplicate public run paths")
    groups = Dict{String,Int}()
    for run in runs
        group = run["group"]
        groups[group] = get(groups, group, 0) + 1
    end
    groups == EXPECTED_GROUPS || error("Public run roster changed: $groups")
    results = [audit_run(root, manifest, run) for run in runs]
    evidence = Dict{String,Any}(
        "schema_version" => 1,
        "manifest" => manifest_path,
        "raw_import_verified_runs" => length(results),
        "kinetic_analysis_ready_runs" => 0,
        "verified_numeric_columns" =>
            sum(length(r["numeric_column_sha256"]) for r in results),
        "numeric_digest_encoding" => "IEEE 754 binary64 little-endian, in acquisition order",
        "rows" => sum(r["rows"] for r in results),
        "segment_count" => sum(length(r["segments"]) for r in results),
        "groups" => groups,
        "runs" => results,
    )
    output = "docs/audits/public_thermal_import.toml"
    open(joinpath(root, output), "w") do io
        return TOML.print(io, evidence; sorted=true)
    end
    println(
        "Verified $(length(results)) public TG/DSC files, $(evidence["rows"]) rows and $(evidence["segment_count"]) segments.",
    )
    return println("Wrote $output")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
