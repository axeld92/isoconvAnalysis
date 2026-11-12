"""Generate known-truth TG CSV controls with independent Python RK4.

Requires Python 3.11+ and its standard library. No Julia routines are used.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tomllib


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "public_tg_energy.toml"
R = 0.00831446261815324


def hazard_path(time, rate, subdivisions):
    """Integrate y=-log(1-alpha), so alpha=-expm1(-y)."""
    y = 0.0
    path = [y]
    for start, end in zip(time, time[1:]):
        h = (end - start) / subdivisions
        for j in range(subdivisions):
            t = start + j * h
            k1 = rate(t, y)
            k2 = rate(t + h / 2, y + h * k1 / 2)
            k3 = rate(t + h / 2, y + h * k2 / 2)
            k4 = rate(t + h, y + h * k3)
            y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        path.append(y)
    return [-math.expm1(-value) for value in path]


def generate(config, control, beta, index, destination):
    cold, hot = config["control_temperature_K"]
    duration = (hot - cold) / beta
    curvature = 0.20 * (-1) ** index if control["curved"] else 0.0
    intercept = control["energy_intercept_kJ_per_mol"]
    slope = control["energy_slope_kJ_per_mol"]
    prefactor = config["control_prefactor_min_inv"]
    rows = config["control_rows"]
    initial_mass, terminal_mass = config["control_mass_mg"]

    def temperature(t):
        phase = t / duration
        return cold + (hot - cold) * (phase + curvature * math.sin(math.tau * phase) / math.tau)

    def derivative(t, hazard):
        energy = intercept - slope * math.expm1(-hazard)
        return prefactor * math.exp(-energy / (R * temperature(t)))

    time = [duration * (i / (rows - 1)) ** 1.18 for i in range(rows)]
    coarse = hazard_path(time, derivative, 4)
    fine = hazard_path(time, derivative, 8)
    difference = max(abs(a - b) for a, b in zip(coarse, fine))
    if difference > config["control_refinement_tolerance_alpha"] or fine[-1] <= max(config["alpha_grid"]):
        raise ValueError(f"Control {control['id']} at {beta} K/min fails trajectory checks")

    run_id = f"{control['id']}_{beta:g}"
    path = destination / f"{run_id}.csv"
    with path.open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("time_min", "temperature_K", "mass_mg"))
        writer.writerows((t, temperature(t), initial_mass - (initial_mass - terminal_mass) * a)
                         for t, a in zip(time, fine))
    return {
        "id": run_id,
        "material_id": control["id"],
        "filename": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "rate_K_per_min": beta,
        "curvature": curvature,
        "refinement_max_abs_alpha": difference,
        "last_alpha": fine[-1],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    destination = parser.parse_args().destination
    destination.mkdir(parents=True, exist_ok=True)
    with CONFIG.open("rb") as stream:
        config = tomllib.load(stream)
    runs = [generate(config, control, beta, index, destination)
            for control in config["controls"]
            for index, beta in enumerate(config["control_rates_K_per_min"])]
    with (destination / "manifest.toml").open("w", encoding="ascii", newline="\n") as stream:
        for run in runs:
            stream.write("[[runs]]\n")
            for key, value in run.items():
                stream.write(f"{key} = {json.dumps(value)}\n")
            stream.write("\n")
    print(f"Generated {len(runs)} independently integrated TG controls.")


if __name__ == "__main__":
    main()
