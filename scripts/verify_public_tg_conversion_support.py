"""Recompute the public TG support audit from the pinned files with Python 3.11+."""

import csv
import hashlib
import io
import math
import statistics
import sys
import tomllib
from pathlib import Path
from typing import NamedTuple


ROOT = Path(__file__).resolve().parents[1]
IMPORT_MANIFEST = "config/public_thermal_import.toml"
TG_PROTOCOL = "config/public_tg_conversion.toml"
DEFAULT_REPORT = ROOT / "docs/audits/public_tg_conversion_support.toml"
PP_10 = "picnik_PP_10_B3"


class Crossing(NamedTuple):
    kind: str
    rows: tuple
    lines: tuple
    temperature_range: tuple
    temperature: float | None
    time: float | None


def require(condition, description):
    if not condition:
        raise ValueError(description)


def close(actual, expected, label, abs_tol=1e-10):
    require(math.isclose(actual, expected, rel_tol=2e-12, abs_tol=abs_tol),
            f"{label} differs: report {actual!r}, raw calculation {expected!r}")


def load_toml(path):
    with path.open("rb") as stream:
        return tomllib.load(stream)


def temperature_window(temperatures, limits):
    rows = [i for i, temperature in enumerate(temperatures)
            if limits[0] <= temperature <= limits[1]]
    require(bool(rows) and rows == list(range(rows[0], rows[-1] + 1)),
            f"Expected a contiguous window for {limits}")
    return slice(rows[0], rows[-1] + 1)


def crossing_events(alpha, temperature, time, lines, first_row, target):
    """Locate the level set of piecewise-linear conversion in acquisition order."""
    result = []
    i = 0
    while i < len(alpha):
        if alpha[i] == target:
            j = i
            while j + 1 < len(alpha) and alpha[j + 1] == target:
                j += 1
            if j > i:
                kind, T, t = "plateau", None, None
            else:
                T, t = temperature[i], time[i]
                if i == 0 or i == len(alpha) - 1:
                    kind = "endpoint"
                elif alpha[i - 1] < target < alpha[i + 1]:
                    kind = "upward"
                elif alpha[i - 1] > target > alpha[i + 1]:
                    kind = "downward"
                else:
                    kind = "tangency"
        elif i + 1 < len(alpha) and min(alpha[i], alpha[i + 1]) < target < max(alpha[i], alpha[i + 1]):
            j = i + 1
            fraction = (target - alpha[i]) / (alpha[j] - alpha[i])
            T = temperature[i] + fraction * (temperature[j] - temperature[i])
            t = time[i] + fraction * (time[j] - time[i])
            kind = "upward" if alpha[j] > alpha[i] else "downward"
        else:
            i += 1
            continue
        bracket = temperature[i:j + 1]
        result.append(Crossing(kind, (i + first_row, j + first_row), (lines[i], lines[j]),
                               (min(bracket), max(bracket)), T, t))
        i = j + 1 if alpha[i] == target else i + 1
    return result


def raw_cases(spec, profile, protocol):
    raw = (ROOT / spec["path"]).read_bytes()
    require(len(raw) == spec["bytes"] and hashlib.sha256(raw).hexdigest() == spec["sha256"],
            f"{spec['id']} raw bytes differ from the pinned input")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter=profile["delimiter"])
    time, temperature_C, mass, lines = [], [], [], []
    for row in reader:
        time.append(float(row[profile["time_column"]]))
        temperature_C.append(float(row[profile["temperature_column"]]))
        mass.append(float(row[profile["signals"][0]["column"]]))
        lines.append(reader.line_num)
    require(len(lines) == spec["rows"] and all(map(math.isfinite, time + temperature_C + mass))
            and all(t1 > t0 for t0, t1 in zip(time, time[1:])),
            f"{spec['id']} has unexpected numeric rows or time order")
    windows = protocol["windows"][spec["group"]]
    initial, terminal, late_initial, early_terminal = (
        statistics.median(mass[temperature_window(temperature_C, windows[name])])
        for name in ("initial_C", "terminal_C", "later_initial_C", "earlier_terminal_C")
    )
    crop = temperature_window(temperature_C, windows["reaction_C"])
    full = slice(0, len(lines))
    cases = {
        "record_endpoints": (full, mass[0], mass[-1]),
        "plateau_medians": (full, initial, terminal),
        "reaction_window": (crop, initial, terminal),
        "later_initial": (full, late_initial, terminal),
        "earlier_terminal": (full, initial, early_terminal),
        "window_and_references": (crop, late_initial, early_terminal),
    }
    if spec["group"] == "experimental_tg":
        cases["zero_terminal"] = (full, initial, 0.0)
    temperature_K = [value + 273.15 for value in temperature_C]
    result = {}
    for name, (rows, m0, mf) in cases.items():
        alpha = [(m0 - value) / (m0 - mf) for value in mass[rows]]
        observed_temperature, observed_time, observed_lines = temperature_K[rows], time[rows], lines[rows]
        result[name] = {
            target: crossing_events(alpha, observed_temperature, observed_time, observed_lines,
                                    rows.start + 1, target)
            for target in protocol["alpha_grid"]
        }
    return result


def crossing_status(events):
    if not events:
        return "missing"
    if len(events) != 1:
        return "ambiguous"
    kind = events[0].kind
    return f"unique_{kind}" if kind in ("upward", "downward") else kind


def check_events(record, events, label):
    require(record["status"] == crossing_status(events) and len(record["events"]) == len(events),
            f"{label}: crossing status or count differs")
    for row, event in zip(record["events"], events, strict=True):
        require((row["kind"], tuple(row["segment_row_range"]), tuple(row["source_line_range"])) ==
                (event.kind, event.rows, event.lines), f"{label}: crossing event or bracket differs")
        for actual, expected in zip(row["temperature_range_K"], event.temperature_range, strict=True):
            close(actual, expected, label + " temperature bracket")
        if event.temperature is not None:
            close(row["temperature_K"], event.temperature, label + " crossing temperature")
            close(row["time_min"], event.time, label + " crossing time")


def check_coverage(record, by_run, target, label):
    require(record["alpha"] == target, f"{label}: unexpected target grid")
    upward = [name for name, events in by_run.items() if crossing_status(events) == "unique_upward"]
    common = len(upward) == len(by_run)
    require((record["upward_run_ids"], record["upward_run_count"], record["common_upward"]) ==
            (upward, len(upward), common), f"{label}: common coverage differs")
    if common:
        temperature = [by_run[name][0].temperature for name in upward]
        cold, hot = min(temperature), max(temperature)
        inverse = [1.0 / value for value in temperature]
        mean_inverse = statistics.fmean(inverse)
        for key, expected, tolerance in (
            ("temperature_min_K", cold, 1e-10),
            ("temperature_mean_K", statistics.fmean(temperature), 1e-10),
            ("temperature_max_K", hot, 1e-10),
            ("temperature_span_K", hot - cold, 1e-10),
            ("reciprocal_temperature_span_K_inv", 1 / cold - 1 / hot, 1e-16),
            ("reciprocal_temperature_Sxx_K_inv_sq",
             math.fsum((value - mean_inverse) ** 2 for value in inverse), 1e-20),
        ):
            close(record[key], expected, f"{label}: {key}", tolerance)
        require(record["temperature_min_run_id"] == upward[temperature.index(cold)] and
                record["temperature_max_run_id"] == upward[temperature.index(hot)],
                f"{label}: extreme temperature run identity differs")
    return common


def verify_points(points, runs, case, targets, raw, events_required):
    common = []
    for target, point in zip(targets, points, strict=True):
        by_run = {run: raw[run][case][target] for run in runs}
        label = f"{case} at alpha={target} with {', '.join(runs)}"
        if check_coverage(point, by_run, target, label):
            common.append(target)
        if events_required:
            for crossing, (run, events) in zip(point["crossings"], by_run.items(), strict=True):
                require(crossing["run_id"] == run, f"{label}: crossing run order differs")
                check_events(crossing, events, f"{run} {case} at alpha={target}")
    return common


def shared_grid(targets, case_grids):
    return [target for target in targets if all(target in grid for grid in case_grids)]


def check_secants(row, reference, central_bounds, raw):
    nodes = [central_bounds[0], 0.50, central_bounds[1]]
    events = [raw[row["id"]][reference][node] for node in nodes]
    if any(crossing_status(event) != "unique_upward" for event in events):
        return
    record = row["central_alpha_program"]
    require(record["crossing_secant_alpha_nodes"] == nodes, f"{row['id']}: secant grid differs")
    events = [event[0] for event in events]
    for actual, left, right in zip(record["crossing_secant_K_per_min"], events[:-1], events[1:], strict=True):
        slope = (right.temperature - left.temperature) / (right.time - left.time)
        close(actual, slope, f"{row['id']}: central conversion crossing secant")


def verify(report_path):
    manifest = load_toml(ROOT / IMPORT_MANIFEST)
    protocol = load_toml(ROOT / TG_PROTOCOL)
    report = load_toml(report_path)
    require(report["raw_import_manifest"] == IMPORT_MANIFEST and report["tg_protocol"] == TG_PROTOCOL
            and report["schema_version"] == 1 and manifest["temperature_unit"] == "C"
            and manifest["time_unit"] == "min", "Unexpected evidence or source schema")
    targets, reference = protocol["alpha_grid"], protocol["comparison_case"]
    require((report["alpha_grid"], report["comparison_case"]) == (targets, reference),
            "Unexpected report grid or reference scenario")
    specs = {run["id"]: run for run in manifest["runs"]
             if run["group"] in ("synthetic_tg", "experimental_tg")}
    raw = {name: raw_cases(run, manifest["profiles"][run["profile"]], protocol)
           for name, run in specs.items()}
    expected_materials = {run["material"] for run in specs.values()}
    require(len(report["groups"]) == len(expected_materials) and
            {group["material_id"] for group in report["groups"]} == expected_materials,
            "Material grouping differs")
    assessments = 0
    group_points = 0
    subset_points = 0
    for group in report["groups"]:
        runs = group["run_ids"]
        expected = {name for name, run in specs.items() if run["material"] == group["material_id"]}
        require(len(runs) == len(expected) and set(runs) == expected, f"Run roster differs: {runs}")
        require([row["id"] for row in group["runs"]] == runs and
                all(specs[row["id"]]["sha256"] == row["source_sha256"] for row in group["runs"]),
                "Evidence source identities or hashes differ")
        for row in group["runs"]:
            check_secants(row, reference, report["central_program_alpha_bounds"], raw)
        cases = set(raw[runs[0]])
        require(len(group["cases"]) == len(cases) and {case["id"] for case in group["cases"]} == cases,
                f"{group['material_id']}: scenario roster differs")
        case_grids, omission_grids = {}, {}
        for case in group["cases"]:
            scenario = case["id"]
            common = verify_points(case["points"], runs, scenario, targets, raw, events_required=True)
            assessments += len(runs) * len(targets)
            group_points += len(targets)
            require(case["common_alpha"] == common, f"{scenario}: common grid differs")
            case_grids[scenario] = common
            expected_omissions = runs if scenario == reference else []
            if scenario != reference and PP_10 in runs:
                expected_omissions = [PP_10]
            omissions = [subset["omitted_run_id"] for subset in case["subsets"]]
            require(omissions == expected_omissions, f"{scenario}: omission roster differs")
            for subset in case["subsets"]:
                selected = [run for run in runs if run != subset["omitted_run_id"]]
                require(subset["run_ids"] == selected, f"{scenario}: subset run IDs differ")
                subset_common = verify_points(subset["points"], selected, scenario, targets, raw,
                                              events_required=False)
                subset_points += len(targets)
                require(subset["common_alpha"] == subset_common and
                        subset["added_vs_full_alpha"] == [a for a in subset_common if a not in common],
                        f"{scenario}: subset common or added grid differs")
                omission_grids.setdefault(subset["omitted_run_id"], []).append(subset_common)
        label = group["material_id"]
        require(group["common_across_all_cases_alpha"] == shared_grid(targets, case_grids.values()),
                f"{label}: shared scenario grid differs")
        expected_intersections = {name: shared_grid(targets, grids) for name, grids in omission_grids.items()
                                  if len(grids) == len(cases)}
        intersections = group["omission_intersections"]
        actual_intersections = {row["omitted_run_id"]: row["common_across_all_cases_alpha"]
                                for row in intersections}
        require(actual_intersections == expected_intersections and len(intersections) == len(actual_intersections),
                f"{label}: shared scenario grid differs for omitted programs")
    require(assessments == report["full_group_crossing_evaluations"] == 1695,
            "Crossing assessment count differs")
    require(subset_points == report["subset_case_count"] * len(targets),
            "Subset point count differs")
    require(group_points == report["group_case_count"] * len(targets), "Group point count differs")
    print(f"Verified {assessments:,} raw crossing assessments, {group_points:,} full-group "
          f"coverage points and {subset_points:,} subset coverage points.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise SystemExit("Usage: python scripts/verify_public_tg_conversion_support.py [report.toml]")
    verify(Path(sys.argv[1]) if len(sys.argv) == 2 else DEFAULT_REPORT)
