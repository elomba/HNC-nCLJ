#!/usr/bin/env python3
"""Analyze dense transport replicas and aggregate partial follow-up results."""
import argparse
import fcntl
import json
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.integrate import cumulative_trapezoid
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent
WINDOWS = [(10, 50), (20, 100), (50, 200), (100, 500), (200, 1000)]
CUTS = [2, 5, 10, 20]
STATES = {
    "B3": "baseline_3C", "U4": "weak_f030_e025",
    "H5": "dispersion_split_d040", "U6": "mono_weak_uniform6C",
    "H7": "mono_weak_split7C",
}
CONTRASTS = {"H5/U4": ("H5", "U4"), "H7/U6": ("H7", "U6"),
             "H5/B3": ("H5", "B3"), "U4/B3": ("U4", "B3")}


def atomic_json(path, value):
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def estimate(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return {"n": 0, "mean": None, "CI": None}
    mean = float(values.mean())
    if len(values) < 2:
        return {"n": 1, "mean": mean, "CI": None}
    error = student_t.ppf(0.975, len(values) - 1) * values.std(ddof=1) / np.sqrt(len(values))
    return {"n": len(values), "mean": mean, "CI": [float(mean - error), float(mean + error)]}


def ratio_estimate(log_ratios):
    result = estimate(log_ratios)
    result["ratio_minus_one"] = None
    result["ratio_CI"] = None
    if result["mean"] is not None:
        result["ratio_minus_one"] = float(np.expm1(result["mean"]))
        result["ratio_CI"] = [float(np.expm1(x)) for x in result["CI"]] if result["CI"] else None
    return result


def welch_ratio(a, b):
    a, b = np.log(np.asarray(a, float)), np.log(np.asarray(b, float))
    if not len(a) or not len(b):
        return {"n_a": len(a), "n_b": len(b), "ratio_minus_one": None, "ratio_CI": None}
    difference = float(a.mean() - b.mean())
    result = {"n_a": len(a), "n_b": len(b), "ratio_minus_one": float(np.expm1(difference)),
              "ratio_CI": None}
    if len(a) > 1 and len(b) > 1:
        va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        se = np.sqrt(va + vb)
        df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
        width = student_t.ppf(0.975, df) * se
        result["ratio_CI"] = [float(np.expm1(difference - width)),
                              float(np.expm1(difference + width))]
    return result


def configurational_z(rows, n_atoms, rho=0.7, temperature=2.0):
    return 1 + (rows[..., 2] - rho * (1 - 1 / n_atoms) * rows[..., 1]) / (
        rho * temperature
    )


def read_positions(folder):
    with Dataset(folder / "positions.nc") as dataset:
        dataset.set_auto_mask(False)
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        typ = np.asarray(dataset.variables["type"][0], dtype=np.int16)
        box = np.asarray(dataset.variables["cell_lengths"][0], dtype=float)
        position = np.asarray(dataset.variables["unwrapped_coordinates"][:], dtype=np.float32)
    assert np.all(np.diff(time) > 0)
    position -= position.mean(axis=1, keepdims=True)
    return time, position, typ, box


def msd_curve(time, position, typ, n_species, max_lag=None, origins=400):
    n_frames = len(time)
    max_lag = min(n_frames // 4, max_lag or n_frames // 4)
    lags = np.unique(np.rint(np.geomspace(1, max_lag, 100)).astype(int))
    time_lags = np.array([np.mean(time[lag:] - time[:-lag]) for lag in lags])
    values = np.empty((len(lags), n_species + 1))
    masks = [typ == i + 1 for i in range(n_species)]
    for j, lag in enumerate(lags):
        start = np.linspace(0, n_frames - lag - 1, min(origins, n_frames - lag),
                            dtype=int)
        squared = np.zeros(len(typ), dtype=np.float64)
        for first in np.array_split(start, max(1, (len(start) + 31) // 32)):
            delta = position[first + lag] - position[first]
            squared += np.sum(np.sum(delta * delta, axis=2), axis=0, dtype=np.float64)
        squared /= len(start)
        values[j, 0] = squared.mean()
        values[j, 1:] = [squared[mask].mean() for mask in masks]
    return time_lags, values


def fit_msd(time_lags, values, windows=WINDOWS):
    fits = {}
    for lo, hi in windows:
        mask = (time_lags >= lo) & (time_lags <= hi)
        if mask.sum() < 5:
            fits[f"{lo}_{hi}"] = None
        else:
            fits[f"{lo}_{hi}"] = [
                float(value) for value in np.polyfit(time_lags[mask], values[mask], 1)[0] / 6
            ]
    return fits


def block_diffusion(time, position, typ, labels, blocks=5):
    results = []
    edges = np.linspace(0, len(time) - 1, blocks + 1, dtype=int)
    for i in range(blocks):
        selected = slice(edges[i], edges[i + 1] + 1)
        lag, curve = msd_curve(time[selected], position[selected], typ, len(labels),
                               origins=100)
        fit = fit_msd(lag, curve, [(20, 100)])["20_100"]
        results.append(dict(zip(["net", *labels], fit)) if fit else None)
    keys = ["net", *labels]
    return {"blocks": results, "summary": {
        key: {"mean": float(np.mean([r[key] for r in results])),
              "standard_error": float(np.std([r[key] for r in results], ddof=1) / np.sqrt(blocks))}
        for key in keys
    }}


def autocorrelation_fft(series, max_lag):
    series = np.asarray(series, dtype=np.float64)
    series -= series.mean(axis=0, keepdims=True)
    n = len(series)
    fft_size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(series, n=fft_size, axis=0)
    correlation = np.fft.irfft(spectrum * spectrum.conjugate(), n=fft_size, axis=0)
    return correlation[:max_lag + 1] / (n - np.arange(max_lag + 1))[:, None]


def viscosity(stress, dt, volume, temperature=2.0):
    max_lag = min(int(round(20 / dt)), len(stress) // 4)
    correlation = autocorrelation_fft(stress, max_lag).mean(axis=1)
    running = cumulative_trapezoid(correlation, dx=dt, initial=0) * volume / temperature
    times = np.arange(max_lag + 1) * dt
    cuts = {str(cut): float(np.interp(cut, times, running)) for cut in CUTS
            if cut <= times[-1]}
    plateau = float(np.mean(running[(times >= 5) & (times <= 10)]))
    return {"cutoffs": cuts, "plateau": plateau}


def yeh_hummer(diffusion, eta, n_atoms, rho=0.7, temperature=2.0):
    if eta <= 0:
        return None
    length = (n_atoms / rho) ** (1 / 3)
    return float(diffusion + 2.837297 * temperature / (6 * np.pi * eta * length))


def vacf(folder, species):
    with Dataset(folder / "velocities.nc") as dataset:
        dataset.set_auto_mask(False)
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        typ = np.asarray(dataset.variables["type"][0], dtype=int)
        velocity = dataset.variables["velocities"]
        com = np.atleast_2d(np.loadtxt(folder / "global_vcom.dat"))
        steps = np.rint(time / 0.002).astype(int)
        _, frame_idx, com_idx = np.intersect1d(steps, com[:, 0].astype(int),
                                                return_indices=True)
        assert len(frame_idx) >= len(time) - 2
        dt = float(np.median(np.diff(time[frame_idx])))
        max_lag = min(int(round(20 / dt)), len(frame_idx) // 4)
        n = len(frame_idx)
        size = 1 << (2 * n - 1).bit_length()
        results = {}
        for index, item in enumerate(species, 1):
            if item["q"] == 0:
                continue
            ids = np.flatnonzero(typ == index)
            assert len(ids)
            curve = np.zeros(max_lag + 1)
            for chunk in np.array_split(ids, max(1, (len(ids) + 31) // 32)):
                velocities = np.asarray(velocity[frame_idx, chunk, :], dtype=np.float32)
                velocities -= com[com_idx, 1:4, None].transpose(0, 2, 1)
                spectrum = np.fft.rfft(velocities, n=size, axis=0)
                curve += np.fft.irfft(spectrum * spectrum.conjugate(), n=size, axis=0)[
                    :max_lag + 1
                ].sum(axis=(1, 2))
            curve /= 3 * len(ids) * (n - np.arange(max_lag + 1))
            integral = cumulative_trapezoid(curve, dx=dt, initial=0)
            results[item["label"]] = {
                str(cut): float(np.interp(cut, np.arange(max_lag + 1) * dt, integral))
                for cut in CUTS if cut <= max_lag * dt
            }
    return results


def concentration_fluctuations(position, typ, box, species):
    frames = np.linspace(0, len(position) - 1, min(len(position), 200), dtype=int)
    result = {}
    for i, item in enumerate(species, 1):
        accum = 0.0
        for frame in frames:
            phase = np.exp(2j * np.pi * position[frame] / box)
            signal = phase[typ == i].sum(axis=0) - item["x"] * phase.sum(axis=0)
            accum += np.mean(np.abs(signal) ** 2)
        result[item["label"]] = float(accum / (len(frames) * len(typ) * item["x"]
                                               * (1 - item["x"])))
    return result


def analyze_replica(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    assert not manifest["smoke"] and (folder / "verified.json").exists()
    species = manifest["composition"]["species"]
    labels = [item["label"] for item in species]
    n_atoms = manifest["N"]
    thermo = np.atleast_2d(np.loadtxt(folder / "thermo_avg.dat"))
    second = thermo[len(thermo) // 2:]
    means = second.mean(axis=0)
    blocks = [chunk.mean(axis=0) for chunk in np.array_split(thermo, 5)]
    z = float(configurational_z(means, n_atoms))
    u = float(means[4] / 2)
    time, position, typ, box = read_positions(folder)
    lag, curve = msd_curve(time, position, typ, len(species))
    fits = fit_msd(lag, curve)
    assert fits["20_100"] is not None
    msd_blocks = block_diffusion(time, position, typ, labels)
    fluctuation = concentration_fluctuations(position, typ, box, species)
    del position
    stress = np.atleast_2d(np.loadtxt(folder / "stress_offdiag.dat"))
    eta = viscosity(stress[:, 1:6], 5 * 0.002, n_atoms / 0.7)
    del stress
    diffusion = dict(zip(["net", *labels], fits["20_100"]))
    diffusion["probe"] = (diffusion["probe_plus2"] + diffusion["probe_minus2"]) / 2
    corrected = {key: yeh_hummer(value, eta["plateau"], n_atoms)
                 for key, value in diffusion.items()}
    result = {
        "name": manifest["composition"]["name"], "campaign": manifest["campaign"],
        "replica": manifest["replica"], "seed": manifest["seed"], "N": n_atoms,
        "thermo": {"Z": z, "u": u, "second_half_rows": len(second),
                   "five_blocks": [{"Z": float(configurational_z(row, n_atoms)),
                                    "u": float(row[4] / 2)} for row in blocks]},
        "transport": {
            "D_net": diffusion["net"], "D_probe": diffusion["probe"],
            "D_species": {label: diffusion[label] for label in labels},
            "D_sign_difference": diffusion["probe_plus2"] - diffusion["probe_minus2"],
            "fitting_windows": fits, "labels": ["net", *labels],
            "block_diffusion": msd_blocks, "Green_Kubo": vacf(folder, species),
            "eta": eta, "D_YH": corrected, "concentration_fluctuations": fluctuation,
            "n_position_frames": len(time),
        },
    }
    np.savetxt(folder / "msd_multi_origin.csv", np.column_stack([lag, curve]), delimiter=",",
               header="lag_tau,net," + ",".join(labels))
    atomic_json(folder / "replica_result.json", result)
    print("Analyzed", folder, flush=True)


def old_rows(name):
    return [json.loads(path.read_text()) for path in sorted(
        (OLD / "md" / name).glob("replica_*/replica_result.json")
    )]


def old_d(row, key):
    transport = row["transport"]
    if key == "net":
        return transport["D_net"]
    if key == "probe":
        return transport["D_probe"]
    return transport["D_species"][key]


def metric(row, key, window="20_100", corrected=False):
    transport = row["transport"]
    if corrected:
        return transport["D_YH"][key]
    fits = transport["fitting_windows"][window]
    if fits is None:
        return None
    labels = transport["labels"]
    if key == "probe":
        return (fits[labels.index("probe_plus2")] + fits[labels.index("probe_minus2")]) / 2
    return fits[labels.index(key)]


def paired_contrast(a, b, key, window="20_100"):
    pairs = sorted(set(a) & set(b))
    values = []
    for index in pairs:
        left, right = metric(a[index], key, window), metric(b[index], key, window)
        if left is not None and right is not None and left > 0 and right > 0:
            values.append(np.log(left / right))
    return ratio_estimate(values)


def aggregate():
    output_dir = ROOT / "results"
    output_dir.mkdir(exist_ok=True)
    with (output_dir / ".aggregate.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rows = [json.loads(path.read_text()) for path in sorted(
            (ROOT / "md").glob("*/*/replica_*/replica_result.json")
        )]
        by = {}
        for row in rows:
            by.setdefault(row["campaign"], {}).setdefault(row["name"], {})[row["replica"]] = row
        systems, contrasts, symmetry = {}, {}, {}
        for campaign in ["long_N4000", "N8000"]:
            systems[campaign], contrasts[campaign], symmetry[campaign] = {}, {}, {}
            group = by.get(campaign, {})
            for code, name in STATES.items():
                reps = group.get(name, {})
                systems[campaign][code] = {
                    key: estimate([value for row in reps.values()
                                   if (value := getter(row)) is not None])
                    for key, getter in {
                        "D_net": lambda r: metric(r, "net"),
                        "D_probe": lambda r: metric(r, "probe"),
                        "D_plus2": lambda r: metric(r, "probe_plus2"),
                        "D_minus2": lambda r: metric(r, "probe_minus2"),
                        "D_sign_difference": lambda r: r["transport"]["D_sign_difference"],
                        "eta": lambda r: r["transport"]["eta"]["plateau"],
                        "D_probe_YH": lambda r: metric(r, "probe", corrected=True),
                        "Z": lambda r: r["thermo"]["Z"],
                        "u": lambda r: r["thermo"]["u"],
                    }.items()
                }
                symmetry[campaign][code] = {}
                for lo, hi in WINDOWS:
                    window = f"{lo}_{hi}"
                    symmetry[campaign][code][window] = estimate([
                        metric(row, "probe_plus2", window) -
                        metric(row, "probe_minus2", window)
                        for row in reps.values()
                        if metric(row, "probe_plus2", window) is not None
                    ])
            for label, (left, right) in CONTRASTS.items():
                contrasts[campaign][label] = {
                    f"{lo}_{hi}": {
                        key: paired_contrast(group.get(STATES[left], {}),
                                             group.get(STATES[right], {}), key, f"{lo}_{hi}")
                        for key in ["net", "probe", "probe_plus2", "probe_minus2"]
                    } for lo, hi in WINDOWS
                }
        finite = {}
        for code, name in STATES.items():
            small = list(by.get("long_N4000", {}).get(name, {}).values())
            large = list(by.get("N8000", {}).get(name, {}).values())
            finite[code] = {}
            for key in ["net", "probe", "probe_plus2", "probe_minus2"]:
                finite[code][key] = {
                    "raw": welch_ratio([metric(row, key) for row in large],
                                       [metric(row, key) for row in small]),
                    "YH": welch_ratio(
                        [metric(row, key, corrected=True) for row in large
                         if metric(row, key, corrected=True) is not None],
                        [metric(row, key, corrected=True) for row in small
                         if metric(row, key, corrected=True) is not None]),
                }
        agreement = {}
        for label in ["H5/U4", "H7/U6"]:
            left = contrasts["long_N4000"][label]["20_100"]["probe"]
            right = contrasts["N8000"][label]["20_100"]["probe"]
            agreement[label] = {
                "same_direction": (None if left["ratio_minus_one"] is None
                                   or right["ratio_minus_one"] is None else
                                   bool(left["ratio_minus_one"] * right["ratio_minus_one"] > 0)),
                "N4000": left, "N8000": right,
            }
        old = {code: {row["replica"]: row for row in old_rows(name)}
               for code, name in STATES.items()}
        combined = {}
        for label, (left, right) in CONTRASTS.items():
            combined[label] = {}
            for key in ["net", "probe"]:
                old_pairs = sorted(set(old[left]) & set(old[right]))
                values = [np.log(old_d(old[left][i], key) / old_d(old[right][i], key))
                          for i in old_pairs]
                new_a = by.get("long_N4000", {}).get(STATES[left], {})
                new_b = by.get("long_N4000", {}).get(STATES[right], {})
                values += [np.log(metric(new_a[i], key) / metric(new_b[i], key))
                           for i in sorted(set(new_a) & set(new_b))]
                combined[label][key] = ratio_estimate(values)
        old_b3 = old["B3"]
        new_b3 = by.get("long_N4000", {}).get(STATES["B3"], {})
        pooled = estimate(
            [old_d(row, "probe_plus2") - old_d(row, "probe_minus2")
             for row in old_b3.values()] +
            [metric(row, "probe_plus2") - metric(row, "probe_minus2")
             for row in new_b3.values()]
        )
        summary = {
            "completed": {"long_N4000": sum(map(len, by.get("long_N4000", {}).values())),
                          "N8000": sum(map(len, by.get("N8000", {}).values()))},
            "expected": {"long_N4000": 30, "N8000": 20},
            "systems": systems, "paired_contrasts": contrasts, "sign_symmetry": symmetry,
            "pooled_B3_sign_difference": {
                "estimate": pooled,
                "note": "Original 3 replicas pooled with new longer 6 replicas; unequal-length runs"
            },
            "finite_size": finite, "contrast_size_agreement": agreement,
            "combined_N4000_secondary": {
                "contrasts": combined,
                "note": "Original 3 plus new 6 independent paired replicas; unequal-length runs"
            },
        }
        atomic_json(output_dir / "summary.json", summary)
        lines = [
            "# Follow-up dense transport",
            "",
            f'Completed: {summary["completed"]["long_N4000"]}/30 long N=4000 and '
            f'{summary["completed"]["N8000"]}/20 N=8000 replicas.',
            "",
            "Primary endpoints: paired H5/U4 and H7/U6 probe diffusion ratios at 20–100 tau; "
            "window stability, B3 charge-sign difference, and N dependence.",
            "",
            "| Campaign | State | n | D_net | D_probe | eta |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for campaign, states in systems.items():
            for code, values in states.items():
                fmt = lambda key, row=values: ("—" if row[key]["mean"] is None
                                               else f'{row[key]["mean"]:.5g}')
                lines.append(f'| {campaign} | {code} | {values["D_net"]["n"]} | '
                             f'{fmt("D_net")} | {fmt("D_probe")} | {fmt("eta")} |')
        lines += ["", "## Primary paired contrasts", "",
                  "| Campaign | Contrast | Probe ratio − 1 | 95% CI |",
                  "|---|---|---:|---:|"]
        for campaign in contrasts:
            for label in ["H5/U4", "H7/U6"]:
                value = contrasts[campaign][label]["20_100"]["probe"]
                point = value.get("ratio_minus_one")
                ci = value.get("ratio_CI")
                lines.append(f'| {campaign} | {label} | '
                             f'{"—" if point is None else f"{point:.3%}"} | '
                             f'{"—" if ci is None else f"[{ci[0]:.3%}, {ci[1]:.3%}]"} |')
        lines += ["", "All estimates and all fitting windows are in summary.json. "
                  "Confidence intervals use independent replicas. The pooled B3 and "
                  "combined N=4000 analyses mix unequal-length runs and are secondary.", ""]
        temp = output_dir / f".report.{os.getpid()}.tmp"
        temp.write_text("\n".join(lines))
        temp.replace(output_dir / "report.md")
        print(f"Aggregated {len(rows)}/50 replicas", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--replica", type=Path)
    args = parser.parse_args()
    if args.replica:
        analyze_replica(args.replica)
    aggregate()
