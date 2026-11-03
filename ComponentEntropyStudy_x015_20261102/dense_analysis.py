#!/usr/bin/env python3
"""Analyze dense transport replicas and aggregate partial follow-up results."""
import json
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.integrate import cumulative_trapezoid

ROOT = Path(__file__).resolve().parent
WINDOWS = [(10, 50), (20, 100), (50, 200), (100, 500), (200, 1000)]
CUTS = [2, 5, 10, 20]


def atomic_json(path, value):
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


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
    region = running[(times >= 5) & (times <= 10)]
    plateau = float(region.mean()) if len(region) else None
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


