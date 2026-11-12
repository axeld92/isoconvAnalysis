"""Independently verify the public PP energy audit with NumPy and SciPy."""

from collections import namedtuple
import csv
import hashlib
import io
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from verify_public_tg_conversion_support import (
    IMPORT_MANIFEST, ROOT, crossing_events, crossing_status, load_toml,
)


R = 0.00831446261815324
NODES, WEIGHTS = np.polynomial.legendre.leggauss(8)
Curve = namedtuple("Curve", "alpha temperature time lines first_row")


def require_close(actual, expected, label, tolerance):
    if not np.allclose(actual, expected, rtol=0, atol=tolerance):
        raise ValueError(f"{label}: expected {expected!r}, got {actual!r}")


def read_pp_runs(report):
    manifest = load_toml(ROOT / IMPORT_MANIFEST)
    specs = {row["id"]: row for row in manifest["runs"]}
    previous = load_toml(ROOT / "docs/audits/public_tg_conversion.toml")
    references = {row["id"]: row["cases"] for row in previous["runs"]}
    runs = {}
    for source in report["polypropylene_runs"]:
        run_id = source["id"]
        spec = specs[run_id]
        profile = manifest["profiles"][spec["profile"]]
        raw = (ROOT / spec["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != spec["sha256"] or source["sha256"] != spec["sha256"]:
            raise ValueError(f"{run_id}: raw source hash mismatch")
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter=profile["delimiter"])
        rows, lines = [], []
        for row in reader:
            rows.append((float(row[profile["time_column"]]), float(row[profile["temperature_column"]]) + 273.15,
                         float(row[profile["signals"][0]["column"]])))
            lines.append(reader.line_num)
        time, temperature, mass = np.array(rows).T
        by_line = {line: index for index, line in enumerate(lines)}
        beta = np.polyfit(time, temperature, 1)[0]
        require_close(beta, source["full_record_OLS_K_per_min"], run_id + " heating rate", 1e-9)
        cases = {}
        for case, metadata in references[run_id].items():
            def mass_reference(reference):
                if reference["method"] == "fixed":
                    return reference["mass_mg"]
                lower, upper = [by_line[line] for line in reference["source_lines"]]
                value = np.median(mass[lower:upper + 1])
                require_close(value, reference["mass_mg"], run_id + " reference mass", 1e-12)
                return value

            initial = mass_reference(metadata["initial_reference"])
            terminal = mass_reference(metadata["terminal_reference"])
            lower, upper = [by_line[line] for line in metadata["reaction_source_lines"]]
            selection = slice(lower, upper + 1)
            alpha = (initial - mass[selection]) / (initial - terminal)
            cases[case] = Curve(alpha, temperature[selection], time[selection],
                                lines[selection], lower + 1)
        runs[run_id] = (beta, cases)
    return runs


def crossing(curve, target):
    events = crossing_events(*curve, target)
    if crossing_status(events) != "unique_upward":
        raise ValueError(f"Expected one crossing at {target}")
    return events[0]


def local_rate(curve, event, width, degree):
    scaled = (curve.time - event.time) / width
    selected = np.abs(scaled) <= 1
    x = scaled[selected]
    weights = np.exp(-0.25 * x * x)
    design = np.vander(x, N=degree + 1, increasing=True)
    coefficients = np.linalg.lstsq(design * weights[:, None], curve.alpha[selected] * weights, rcond=None)[0]
    return coefficients[1] / width


def gauss_path(curve, event, target, delta):
    start = crossing(curve, target - delta)
    inside = (curve.time > start.time) & (curve.time < event.time)
    if not np.all(np.diff(np.r_[target - delta, curve.alpha[inside], target]) > 0):
        raise ValueError(f"Nonmonotone conversion in the interval ending at {target}")
    t = np.r_[start.time, curve.time[inside], event.time]
    T = np.r_[start.temperature, curve.temperature[inside], event.temperature]
    temperatures = (T[:-1] + T[1:])[:, None] / 2 + np.diff(T)[:, None] * NODES / 2
    weights = np.diff(t)[:, None] * WEIGHTS / 2
    return temperatures, weights


def nonlinear_energy(curves, events, target, settings):
    paths = [gauss_path(curve, event, target, settings["delta_alpha"])
             for curve, event in zip(curves, events, strict=True)]

    def objective(energy):
        logs = np.array([np.log(np.sum(weight * np.exp(-energy / (R * temp))))
                         for temp, weight in paths])
        # Symmetric pair excess is evaluated without subtracting nearly equal numbers.
        return np.sum(2 * np.sinh((logs[:, None] - logs[None, :]) / 2) ** 2)

    fit = minimize_scalar(objective, bounds=settings["energy_bounds_kJ_per_mol"],
                          method="bounded", options={"xatol": 1e-8})
    if not fit.success:
        raise ValueError(f"Independent energy minimization failed at {target}")
    return fit.x


def verify_profile(profile, runs):
    method, case, settings = (profile[name] for name in ("method", "case", "settings"))
    curves = [runs[name][1][case] for name in profile["run_ids"]]
    maximum_error = 0.0
    for point in profile["points"]:
        target = point["alpha"]
        label = f"{method}/{case}/alpha={target}"
        events = [crossing(curve, target) for curve in curves]
        T = np.array([event.temperature for event in events])
        require_close(T, point["crossing_temperature_K"], label + " crossing temperatures", 1e-9)
        if [list(event.lines) for event in events] != point["crossing_source_lines"]:
            raise ValueError(f"{label}: crossing provenance differs")
        if method == "advanced_vyazovkin":
            estimate = nonlinear_energy(curves, events, target, settings)
            tolerance = 0.005
        else:
            if method == "friedman":
                rates = [local_rate(curve, event, width, settings["degree"])
                         for curve, event, width in zip(curves, events, settings["half_windows_min"], strict=True)]
                power, coefficient = 0.0, 1.0
            else:
                rates = [runs[name][0] for name in profile["run_ids"]]
                power, coefficient = {"kas": (2.0, 1.0), "fwo": (0.0, 1.052),
                                      "starink": (1.92, 1.0008)}[method]
            require_close(rates, point["method_inputs"], label + " rates", 1e-9)
            estimate = -R * np.polyfit(1 / T, np.log(rates / T ** power), 1)[0] / coefficient
            tolerance = 1e-6
        require_close(estimate, point["energy_kJ_per_mol"], label + " energy", tolerance)
        maximum_error = max(maximum_error, abs(estimate - point["energy_kJ_per_mol"]))
    return len(profile["points"]), maximum_error


def main():
    report = load_toml(ROOT / "docs/audits/public_tg_energy.toml")
    for path, fingerprint in report["input_sha256"].items():
        if hashlib.sha256((ROOT / Path(path)).read_bytes()).hexdigest() != fingerprint:
            raise ValueError(f"Changed audit input {path}")
    runs = read_pp_runs(report)
    total = 0
    maximum_error = 0.0
    for cohort in report["polypropylene_cohorts"]:
        for profile in [*cohort["profiles"], cohort["refined_primary"]]:
            count, difference = verify_profile(profile, runs)
            total += count
            maximum_error = max(maximum_error, difference)
    print(f"Verified {total} PP estimates. Maximum independent difference {maximum_error:.6g} kJ/mol.")


if __name__ == "__main__":
    main()
