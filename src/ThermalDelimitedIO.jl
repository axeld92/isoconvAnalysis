const _CP1252_EXTENDED = (
    '€',
    '\0',
    '‚',
    'ƒ',
    '„',
    '…',
    '†',
    '‡',
    'ˆ',
    '‰',
    'Š',
    '‹',
    'Œ',
    '\0',
    'Ž',
    '\0',
    '\0',
    '‘',
    '’',
    '“',
    '”',
    '•',
    '–',
    '—',
    '˜',
    '™',
    'š',
    '›',
    'œ',
    '\0',
    'ž',
    'Ÿ',
)

function _decode_thermal_export(bytes, encoding, fail)
    if encoding === :utf8
        text = String(bytes)
        isvalid(text) || fail("invalid UTF-8. Supply the export encoding explicitly")
        return startswith(text, '\ufeff') ? chop(text; head=1, tail=0) : text
    end
    return sprint() do io
        for (offset, byte) in enumerate(bytes)
            char = 0x80 <= byte <= 0x9f ? _CP1252_EXTENDED[Int(byte) - 127] : Char(byte)
            char == '\0' && fail("invalid Windows-1252 text at byte $offset")
            print(io, char)
        end
    end
end

function _thermal_delimited_fields(line, delimiter, fail)
    occursin('"', line) || return strip.(split(line, delimiter; keepempty=true))
    fields = String[]
    field = IOBuffer()
    state = :start
    for char in line
        if state === :quoted
            if char == '"'
                state = :after_quote
            else
                print(field, char)
            end
        elseif state === :after_quote && char == '"'
            print(field, char)
            state = :quoted
        elseif char == delimiter
            value = String(take!(field))
            push!(fields, state in (:after_quote, :after_space) ? value : strip(value))
            state = :start
        elseif state in (:after_quote, :after_space)
            isspace(char) || fail("unexpected text after a closing quote")
            state = :after_space
        elseif char == '"'
            state === :start || fail("unexpected quote in an unquoted field")
            state = :quoted
        elseif state !== :start || !isspace(char)
            print(field, char)
            state = :bare
        end
    end
    state === :quoted && fail("unterminated quoted field. Multiline fields are unsupported")
    value = String(take!(field))
    push!(fields, state in (:after_quote, :after_space) ? value : strip(value))
    return fields
end

function _thermal_metadata!(metadata, line, delimiter, fail)
    startswith(line, '#') || fail("expected a #key$(delimiter)value record or the header")
    fields = _thermal_delimited_fields(chop(line; head=1, tail=0), delimiter, fail)
    length(fields) == 2 || fail("expected exactly two metadata fields")
    key = strip(endswith(fields[1], ':') ? chop(fields[1]; tail=1) : fields[1])
    isempty(key) && fail("empty metadata key")
    haskey(metadata, key) && fail("duplicate metadata key $(repr(key))")
    metadata[key] = fields[2]
    return nothing
end

function _thermal_numeric(fields, index, name, fail)
    value = tryparse(Float64, fields[index])
    if isnothing(value) || !isfinite(value)
        fail("expected a finite number in column $(repr(name)), got $(repr(fields[index]))")
    end
    return value
end

"""
    read_thermal_delimited(path, mapping; id, material_id, delimiter=',',
                           encoding=:utf8, header_prefix=nothing, unit_labels=nothing,
                           expected_sha256=nothing)

Import thermal observations from a delimited file with an explicit `ThermalColumnMapping`.
`id` defaults to the filename stem. With `header_prefix`, the first nonblank line starting
with that prefix gives the column names and preceding `#key<delimiter>value` records give
metadata. Without a prefix, the first nonblank line is the header. The prefix is removed
before parsing header fields. Column names are trimmed. Records may use LF or CRLF line
endings, single-line quoted fields, and doubled quotes.
Choose `encoding=:cp1252` explicitly for Windows-1252 exports. An initial UTF-8 BOM is
accepted with `encoding=:utf8`. `expected_sha256` optionally verifies the exact source bytes.
Pass `unit_labels` to check a row of literal unit labels after the column names. The
labels follow the file's column order and the mapping still specifies the units used.

Each data row must have the header's field count and finite values in all mapped numeric
columns. A mapped segment column supplies nonempty labels in source order. Unmapped
columns are recorded by name in the provenance. Units, signal sign conventions, and
segment labels are taken from the explicit mapping and source columns.
"""
function read_thermal_delimited(
    path::AbstractString,
    mapping::ThermalColumnMapping;
    id=splitext(basename(path))[1],
    material_id=nothing,
    delimiter::Char=(','),
    encoding::Symbol=:utf8,
    header_prefix::Union{Nothing,AbstractString}=nothing,
    unit_labels::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    expected_sha256::Union{Nothing,AbstractString}=nothing,
)
    fail(message) = throw(DataImportError(String(path), "", message))
    if delimiter in ('"', '\0') || (isspace(delimiter) && delimiter != '\t')
        fail("invalid delimiter $(repr(delimiter))")
    end
    encoding in (:utf8, :cp1252) || fail("unsupported encoding $(repr(encoding))")
    if !isnothing(header_prefix) && isempty(header_prefix)
        fail("header_prefix must be nonempty")
    end
    bytes = try
        read(path)
    catch err
        err isa Union{SystemError,Base.IOError} || rethrow()
        fail("could not read file: $(sprint(showerror, err))")
    end
    sha256 = bytes2hex(SHA.sha256(bytes))
    if !isnothing(expected_sha256)
        occursin(r"\A[0-9a-fA-F]{64}\z", expected_sha256) ||
            fail("invalid expected SHA-256")
        sha256 == lowercase(expected_sha256) ||
            fail("SHA-256 mismatch: expected $expected_sha256, got $sha256")
    end
    source = _decode_thermal_export(bytes, encoding, fail)
    occursin('\0', source) && fail("source contains a NUL character")
    return _parse_thermal_export(
        source,
        path,
        mapping,
        id,
        material_id,
        delimiter,
        header_prefix,
        unit_labels,
        sha256,
        fail,
    )
end

function _parse_thermal_export(
    source,
    path,
    mapping,
    id,
    material_id,
    delimiter,
    header_prefix,
    unit_labels,
    sha256,
    fail,
)
    metadata = Dict{String,String}()
    source_lines = Int[]
    raw_temperature = Float64[]
    raw_time = isnothing(mapping.time) ? nothing : Float64[]
    raw_signals = [Float64[] for signal in mapping.signals]
    labels = isnothing(mapping.segment_column) ? nothing : String[]
    temperature_name = mapping.temperature.name
    time_name = isnothing(mapping.time) ? nothing : mapping.time.name
    signal_names = [signal.column.name for signal in mapping.signals]
    segment_name = mapping.segment_column
    header = String[]
    unmapped = String[]
    temperature_index = time_index = segment_index = 0
    signal_indices = Int[]
    units_pending = !isnothing(unit_labels)
    for (line_number, line) in enumerate(eachline(IOBuffer(source)))
        isempty(strip(line)) && !occursin(delimiter, line) && continue
        fail_line(message) = fail("line $line_number: $message")
        if isempty(header)
            if !isnothing(header_prefix) && !startswith(line, header_prefix)
                _thermal_metadata!(metadata, line, delimiter, fail_line)
                continue
            end
            header_line = if isnothing(header_prefix)
                line
            else
                chop(line; head=length(header_prefix), tail=0)
            end
            fields = _thermal_delimited_fields(header_line, delimiter, fail_line)
            header = String.(strip.(fields))
            any(isempty, header) && fail_line("empty column name")
            allunique(header) || fail_line("duplicate column names")
            index = Dict(name => position for (position, name) in enumerate(header))
            require_column(name) = get(index, name) do
                return fail_line("missing mapped column $(repr(name))")
            end
            temperature_index = require_column(temperature_name)
            signal_indices = require_column.(signal_names)
            isnothing(time_name) || (time_index = require_column(time_name))
            isnothing(segment_name) || (segment_index = require_column(segment_name))
            mapped = Set([temperature_index, time_index, segment_index, signal_indices...])
            unmapped = [name for (position, name) in enumerate(header) if position ∉ mapped]
            continue
        end
        fields = _thermal_delimited_fields(line, delimiter, fail_line)
        length(fields) == length(header) ||
            fail_line("expected $(length(header)) fields, got $(length(fields))")
        if units_pending
            strip.(fields) == unit_labels || fail_line("column unit labels do not match")
            units_pending = false
            continue
        end
        push!(
            raw_temperature,
            _thermal_numeric(fields, temperature_index, temperature_name, fail_line),
        )
        if !isnothing(raw_time)
            push!(raw_time, _thermal_numeric(fields, time_index, time_name, fail_line))
        end
        for (values, index, name) in zip(raw_signals, signal_indices, signal_names)
            push!(values, _thermal_numeric(fields, index, name, fail_line))
        end
        if !isnothing(labels)
            label = fields[segment_index]
            isempty(strip(label)) &&
                fail_line("empty segment label in $(repr(segment_name))")
            push!(labels, label)
        end
        push!(source_lines, line_number)
    end
    isempty(header) && fail("missing header")
    units_pending && fail("missing column unit labels")
    isempty(source_lines) && fail("no data rows")
    provenance = ThermalProvenance(abspath(path), sha256, source_lines, unmapped)
    return _build_thermal_run(
        id,
        material_id,
        mapping,
        raw_temperature,
        raw_time,
        raw_signals,
        labels,
        metadata,
        provenance,
    )
end
