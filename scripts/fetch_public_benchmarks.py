"""Acquire immutable public benchmark samples; never execute upstream source code.

Default source selection is small: mixchar's two examples, pICNIK's synthetic/PP
examples, and the five dynamic plus seven isothermal EU TURBO DSC CSVs. Upstream
bytes are retained without normalization. --offline verifies already saved bytes.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import statistics
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "external"
AGENT = "isoconvAnalysis-public-benchmark-audit/1.0"
REPOSITORIES = {
    "mixchar": ("smwindecker/mixchar", "8e26d64d6603e583e6cbcddb1ea6123cd130963f"),
    "picnik": ("ErickErock/pICNIK", "9156453c65fbc6d540eacb6ff371309eb2e5a11c"),
    "cast_pmma": ("firedynamics/FireSciPy", "71585a5e889a405e8fb3ab9b1d8a2243766e37e4"),
}
PMMA_FOLDER = "docs/tutorials/pyrolysis/data/"
PMMA_RAMP_FILE = re.compile(r"tga_dynamic_n2_dyn(?:5|10|20|30|60)_powder_1mg_r[1-3]\.txt")
PMMA_DOIS = ("10.24355/dbbs.084-202504170956-0", "10.24355/dbbs.084-202508070641-0")


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": AGENT})
    with urllib.request.urlopen(req, timeout=90) as response:
        return response.read()


def json_fetch(url):
    return json.loads(fetch(url))


def digest(payload, algorithm="sha256"):
    return hashlib.new(algorithm, payload).hexdigest()


def store_payload(directory, name, url, expected=None):
    path = directory / "raw" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = path.read_bytes() if path.exists() else fetch(url)
    if expected:
        algorithm, value = expected.split(":", 1)
        if digest(payload, algorithm) != value:
            raise ValueError(f"Checksum mismatch for {url}")
    if not path.exists():
        path.write_bytes(payload)
    return {"path": str(path.relative_to(ROOT)).replace("\\", "/"),
            "url": url, "bytes": len(payload), "sha256": digest(payload),
            "upstream_checksum": expected or "not supplied"}


def acquire_github(source):
    repo, revision = REPOSITORIES[source]
    directory = OUT / source
    directory.mkdir(parents=True, exist_ok=True)
    tree_url = f"https://api.github.com/repos/{repo}/git/trees/{revision}?recursive=1"
    tree = json_fetch(tree_url)
    if tree.get("truncated"):
        raise ValueError("Truncated GitHub inventory")
    metadata = {"repository": repo, "revision": revision, "tree": tree,
                "metadata_url": tree_url}
    files = []
    for entry in tree["tree"]:
        p = entry["path"]
        if source == "mixchar":
            selected = p in ("data/juncus.rda", "data/marsilea.rda", "LICENSE", "DESCRIPTION",
                             "R/juncus-data.R", "R/marsilea-data.R", "vignettes/mixchar.Rmd")
        elif source == "picnik":
            selected = (p.endswith(".csv") and p.startswith(("examples/Constant_E/", "examples/Two_Steps/",
                        "examples/Variable_E/", "examples/PolyProp/")) and ".ipynb_checkpoints" not in p)
            selected |= p in ("LICENSE", "README.md", "examples/example.py")
        else:
            name = p.removeprefix(PMMA_FOLDER)
            selected = p.startswith(PMMA_FOLDER) and (name == "README.txt" or PMMA_RAMP_FILE.fullmatch(name))
        if selected:
            url = f"https://raw.githubusercontent.com/{repo}/{revision}/" + urllib.parse.quote(p)
            record = store_payload(directory, p, url)
            if record["bytes"] != entry["size"]:
                raise ValueError(f"GitHub byte count mismatch: {p}")
            payload = (ROOT / record["path"]).read_bytes()
            git_blob = hashlib.sha1(f"blob {len(payload)}\0".encode() + payload).hexdigest()
            if git_blob != entry["sha"]:
                raise ValueError(f"Git blob checksum mismatch: {p}")
            record["upstream_git_blob_sha1"] = git_blob
            files.append(record)
    if source == "cast_pmma":
        if len(files) != 16:
            raise ValueError("Expected fifteen cast PMMA runs and their README")
        records = {doi: json_fetch(f"https://api.datacite.org/dois/{doi}")["data"]["attributes"]
                   for doi in PMMA_DOIS}
        if any(not any(r.get("rightsIdentifier") == "cc-by-4.0" for r in m["rightsList"])
               for m in records.values()):
            raise ValueError("Expected CC BY 4.0 on both DOI records")
        metadata = {"repository": repo, "revision": revision, "metadata_url": tree_url,
                    "source_doi": PMMA_DOIS[0], "datacite_records": records,
                    "data_license_status": "CC BY 4.0 (original cast PMMA deposition)"}
    else:
        metadata["code_license"] = "MIT + file LICENSE (DESCRIPTION)" if source == "mixchar" else "MIT (repository LICENSE)"
        metadata["data_license_status"] = "No separate dataset-specific license found; local evaluation only pending redistribution review"
    (directory / "source_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return files, metadata


def acquire_turbo():
    directory = OUT / "zenodo_8113838"
    directory.mkdir(parents=True, exist_ok=True)
    url = "https://zenodo.org/api/records/8113838"
    metadata = json_fetch(url)
    if metadata["metadata"]["license"]["id"] != "cc-by-4.0":
        raise ValueError("Expected the reviewed CC-BY-4.0 data license")
    files = []
    for item in metadata["files"]:
        if item["key"].startswith(("Dynamic_", "Isotherm_")):
            record = store_payload(directory, item["key"], item["links"]["self"], item["checksum"])
            if record["bytes"] != item["size"]:
                raise ValueError("Zenodo byte count mismatch")
            files.append(record)
    if len(files) != 12:
        raise ValueError("Expected five dynamic and seven isothermal CSV files")
    (directory / "source_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return files, metadata


def text_csv_inventory(path):
    """Structural inventory only; retains units/sign/conversion decisions for audit."""
    payload = path.read_bytes()
    encodings = ("utf-16",) if payload.startswith((b"\xff\xfe", b"\xfe\xff")) else ("utf-8-sig", "cp1252")
    for encoding in encodings:
        try:
            text = payload.decode(encoding)
            break
        except UnicodeError:
            continue
    else:
        return {"parse_status": "text decoding failed"}
    lines = text.splitlines()
    preview = lines[:12]
    delimiter = max((",", ";", "\t"), key=lambda s: sum(line.count(s) for line in lines[:50]))
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    summary = {"parse_status": "structural inventory only", "encoding": encoding,
            "delimiter": repr(delimiter), "line_count": len(lines), "first_lines": preview,
            "row_width_counts": {str(n): sum(len(r) == n for r in rows) for n in sorted({len(r) for r in rows})}}
    header_index = next((i for i, row in enumerate(rows) if row and row[0].startswith("##Temp.")), None)
    if header_index is None:
        header_index = 0
    header = [v.strip() for v in rows[header_index]]
    numeric, rejected = [], 0
    for row in rows[header_index + 1:]:
        if not row or all(not value.strip() for value in row):
            continue
        try:
            values = [float(value) for value in row]
            if len(values) != len(header) or not all(math.isfinite(value) for value in values):
                raise ValueError("nonfinite or wrong width")
            numeric.append(values)
        except ValueError:
            rejected += 1
    summary.update({"column_names": header, "numeric_rows": len(numeric), "rejected_rows": rejected})
    if not numeric:
        return summary
    temperature_index = next(i for i, name in enumerate(header) if "temp" in name.lower() or name.startswith("T ("))
    time_index = next(i for i, name in enumerate(header) if "time" in name.lower() or name.startswith("t ("))
    signal_index = next(i for i, name in enumerate(header) if "DSC" in name or "Weight" in name or name.startswith("m ("))
    segment_index = next((i for i, name in enumerate(header) if name == "Segment"), None)
    metadata = {row[0].lstrip("#").strip().rstrip(":"): row[1].strip() for row in rows[:header_index] if len(row) == 2}
    summary["instrument_metadata"] = metadata
    summary["segments"] = []
    groups = sorted({int(row[segment_index]) if segment_index is not None else 0 for row in numeric})
    for segment in groups:
        subset = [row for row in numeric if segment_index is None or int(row[segment_index]) == segment]
        times = [row[time_index] for row in subset]
        temperatures = [row[temperature_index] for row in subset]
        signals = [row[signal_index] for row in subset]
        mean_t, mean_T = statistics.mean(times), statistics.mean(temperatures)
        ss = sum((t - mean_t) ** 2 for t in times)
        slope = sum((t - mean_t) * (T - mean_T) for t, T in zip(times, temperatures)) / ss if ss else None
        summary["segments"].append({"segment": segment, "rows": len(subset),
            "time_range_reported": [min(times), max(times)], "temperature_range_reported": [min(temperatures), max(temperatures)],
            "signal_range_reported": [min(signals), max(signals)],
            "nonincreasing_time_steps": sum(b <= a for a, b in zip(times, times[1:])),
            "temperature_decrease_steps": sum(b < a for a, b in zip(temperatures, temperatures[1:])),
            "linear_temperature_time_slope_reported_units": slope})
    return summary


def verify_saved(source):
    directory = OUT / source
    manifest = json.loads((directory / "acquisition_manifest.json").read_text(encoding="utf-8"))
    for record in manifest["files"]:
        path = ROOT / record["path"]
        if digest(path.read_bytes()) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ValueError(f"Changed raw payload: {path}")
    return len(manifest["files"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", nargs="+", choices=[*REPOSITORIES, "zenodo_8113838"],
                        default=["mixchar", "picnik", "zenodo_8113838"])
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    for source in args.sources:
        if args.offline:
            print(source, "verified files:", verify_saved(source))
            continue
        files, metadata = acquire_turbo() if source == "zenodo_8113838" else acquire_github(source)
        directory = OUT / source
        if source != "cast_pmma":
            inventory = {r["path"]: text_csv_inventory(ROOT / r["path"]) for r in files if r["path"].endswith(".csv")}
            (directory / "csv_inventory.json").write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8")
        manifest = {"source_id": source, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    "retrieval_script_sha256": digest(Path(__file__).read_bytes()),
                    "normalization_applied": False, "upstream_code_executed": False,
                    "files": files}
        (directory / "acquisition_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(source, "acquired files:", len(files), "bytes:", sum(r["bytes"] for r in files))


if __name__ == "__main__":
    main()
