using Downloads
using SHA
using TOML

const PUBLIC_TG_ROOT = normpath(joinpath(@__DIR__, ".."))
const PICNIK_RAW_PREFIX = "data/external/picnik/raw/"
const PICNIK_URL =
    "https://raw.githubusercontent.com/ErickErock/pICNIK/" *
    "9156453c65fbc6d540eacb6ff371309eb2e5a11c/"

function picnik_source_url(path)
    startswith(path, PICNIK_RAW_PREFIX) || error("Expected a pICNIK raw path: $path")
    relative = chop(path; head=length(PICNIK_RAW_PREFIX), tail=0)
    parts = split(relative, '/')
    all(part -> occursin(r"^[A-Za-z0-9_. -]+$", part) && part ∉ (".", ".."), parts) ||
        error("Invalid pICNIK path: $path")
    return PICNIK_URL * replace(relative, " " => "%20")
end

function verify_public_tg_file(path, spec)
    filesize(path) == spec["bytes"] || error("Byte count mismatch: $path")
    bytes2hex(open(sha256, path)) == spec["sha256"] || error("SHA-256 mismatch: $path")
    return nothing
end

function fetch_public_tg(root; offline=false)
    manifest = TOML.parsefile(joinpath(root, "config", "public_thermal_import.toml"))
    specs = filter(r -> r["group"] in ("experimental_tg", "synthetic_tg"), manifest["runs"])
    paths = getindex.(specs, "path")
    length(paths) == 18 && allunique(paths) || error("Unexpected public TG run roster")

    for spec in specs
        relative = spec["path"]
        url = picnik_source_url(relative)
        path = joinpath(root, split(relative, '/')...)
        if isfile(path)
            verify_public_tg_file(path, spec)
            continue
        end
        offline && error("Missing cached pICNIK input: $relative")
        mkpath(dirname(path))
        mktemp(dirname(path)) do temporary, io
            close(io)
            Downloads.download(url, temporary; timeout=90)
            verify_public_tg_file(temporary, spec)
            return mv(temporary, path)
        end
    end
    println("Verified $(length(paths)) pinned pICNIK TG files.")
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    ARGS in (String[], ["--offline"]) ||
        error("Usage: julia --project=. scripts/fetch_public_tg.jl [--offline]")
    fetch_public_tg(PUBLIC_TG_ROOT; offline=(!isempty(ARGS)))
end
