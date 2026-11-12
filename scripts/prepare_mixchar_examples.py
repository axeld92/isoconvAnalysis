"""Extract the two pinned mixchar R data.frames to local CSV files.

Run with Python 3.11+ and the optional decoder installed in an isolated Python
environment (`python -m pip install --only-binary=:all: rdata==1.0.0`).
Only byte-verified inputs from config/mixchar_thermal_import.toml are decoded.
"""

from hashlib import sha256
from io import BytesIO, StringIO
from math import isfinite
from pathlib import Path
from tempfile import NamedTemporaryFile
import csv
import struct
import tomllib


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "mixchar_thermal_import.toml"


def verify_hash(data, expected, label):
    actual = sha256(data).hexdigest()
    if actual != expected:
        raise ValueError(f"{label}: expected SHA-256 {expected}, got {actual}")


def extract_table(rdata, dataframe_type, run, columns):
    path = ROOT / run["rda_path"]
    source = path.read_bytes()
    if len(source) != run["rda_bytes"]:
        raise ValueError(f"{path}: expected {run['rda_bytes']} bytes, got {len(source)}")
    verify_hash(source, run["rda_sha256"], path)
    name = run["id"]
    objects = rdata.read_rda(
        BytesIO(source),
        constructor_dict={"data.frame": rdata.conversion.dataframe_constructor},
    )
    if set(objects) != {name} or not isinstance(objects[name], dataframe_type):
        raise ValueError(f"{path}: expected one R data.frame named {name}")
    frame = objects[name]
    if tuple(frame.columns) != columns or frame.shape != (run["rows"], len(columns)):
        raise ValueError(f"{path}: unexpected data.frame columns or dimensions")
    if list(map(str, frame.dtypes)) != ["float64"] * len(columns):
        raise ValueError(f"{path}: expected double-precision numeric columns")
    return frame


def verified_csv(frame, run, columns):
    expected_hashes = (run["temperature_sha256_ieee754_le"], run["signal_sha256_ieee754_le"])
    for column, expected in zip(columns, expected_hashes, strict=True):
        values = frame[column].tolist()
        if not all(map(isfinite, values)):
            raise ValueError(f"{run['id']}: nonfinite data in {column}")
        packed_values = struct.pack(f"<{len(values)}d", *values)
        verify_hash(packed_values, expected, f"{run['id']} {column}")

    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(columns)
    writer.writerows(
        tuple(repr(float(value)) for value in row)
        for row in frame.itertuples(index=False, name=None)
    )
    data = output.getvalue().encode("utf-8")
    verify_hash(data, run["csv_sha256"], run["csv_path"])
    return data


def write_csv(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    output = NamedTemporaryFile(dir=path.parent, prefix=f".{path.stem}-", delete=False)
    temporary_path = Path(output.name)
    try:
        with output:
            output.write(data)
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)


def main():
    with CONFIG.open("rb") as config_file:
        config = tomllib.load(config_file)
    version = config["decoder_version"]
    try:
        import rdata
        from pandas import DataFrame
    except ImportError as error:
        raise SystemExit(
            f"Install rdata=={version} in an isolated Python environment. {error}"
        ) from error
    if rdata.__version__ != version:
        raise SystemExit(f"Expected rdata {version}, found {rdata.__version__}.")

    columns = (config["temperature_column"], config["signal_column"])
    outputs = []
    for run in config["runs"]:
        frame = extract_table(rdata, DataFrame, run, columns)
        outputs.append((run, verified_csv(frame, run, columns)))
    for run, data in outputs:
        write_csv(ROOT / run["csv_path"], data)
        print(f"{run['csv_path']}: verified {run['rows']} rows, SHA-256 {run['csv_sha256']}")


if __name__ == "__main__":
    main()
