#!/usr/bin/env python3
"""Analyze one high-divalent MD allocation, including charge-current transport."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator

from dense_analysis import (WINDOWS, autocorrelation_fft, atomic_json,
                            block_diffusion, concentration_fluctuations,
                            configurational_z, fit_msd, msd_curve,
                            read_positions, vacf, viscosity, yeh_hummer)

ROOT = Path(__file__).resolve().parent
N = 4000
T = 2.0
RHO = 0.7
DT = 0.002


def entropy(rho, z, pe, q2, pchip=False):
    """TI excess entropy with the Debye-Huckel sqrt(rho) endpoint."""
    x = np.sqrt(rho)
    coefficient = (4 * np.pi * 0.5 * q2) ** 1.5 / (12 * np.pi)
    interpolator = PchipInterpolator if pchip else CubicSpline
    integral = interpolator(np.r_[0, x], np.r_[-coefficient, 2 * (z - 1) / x]).antiderivative()
    return pe / T - (integral(x) - integral(0))


def charge_current(charges, velocities):
    """Sum q*v; a neutral system is invariant to a uniform velocity shift."""
    charges = np.asarray(charges, float)
    assert abs(charges.sum()) < 1e-10
    return np.einsum('i,...ij->...j', charges, velocities)


def conductivity_integral(correlation, dt, volume, temperature=T):
    """Integrate the dot-current autocorrelation with 1/(3Vk_BT)."""
    return cumulative_trapezoid(correlation, dx=dt, initial=0) / (3 * volume * temperature)


def conductivity(current, dt, volume, temperature=T):
    """GK conductivity in LJ units, using LAMMPS units lj (Coulomb prefactor 1).

    Charges are LAMMPS dimensionless q, velocity is sigma/tau, and k_B=1.
    sigma has units q^2/(epsilon*tau*sigma). The formula includes 1/(3VT).
    """
    current = np.asarray(current, float)
    assert current.ndim == 2 and current.shape[1] == 3 and len(current) > 4
    max_lag = min(round(20 / dt), len(current) // 4)
    correlation = autocorrelation_fft(current.copy(), max_lag).sum(axis=1)
    running = conductivity_integral(correlation, dt, volume, temperature)
    time = np.arange(max_lag + 1) * dt
    cuts = {str(cut): float(np.interp(cut, time, running))
            for cut in [2, 5, 10, 20] if cut <= time[-1]}
    region = running[(time >= 5) & (time <= 10)]
    return {"cutoffs": cuts, "plateau": float(region.mean()) if len(region) else None,
            "max_lag_tau": float(time[-1]), "time": time, "running": running,
            "correlation": correlation}


def nernst_einstein(species, diffusion, rho=RHO, temperature=T):
    """N/(Vk_BT) sum_i x_i q_i^2 D_i, in the same LJ conductivity units."""
    return float(rho / temperature * sum(
        item['x'] * item['q'] ** 2 * diffusion[item['label']] for item in species
    ))


def thermo_rows(folder, manifest):
    rhos = np.asarray(manifest['densities'], float)
    rows = []
    for rho in rhos:
        stage = folder / ('transport' if np.isclose(rho, RHO) else f'rho_{rho:.4f}')
        assert (stage / 'verified.json').is_file(), stage
        data = np.atleast_2d(np.loadtxt(stage / 'thermo_avg.dat'))
        assert data.shape[1] == 5 and len(data) > 0
        rows.append(data[len(data) // 2:].mean(axis=0))
    return rhos, np.asarray(rows)


def analyze(folder, smoke=False):
    manifest = json.loads((folder / 'manifest.json').read_text())
    assert bool(manifest['smoke']) == smoke
    species = manifest['composition']['species']
    assert sum(round(N * s['x']) * s['q'] for s in species) == 0
    rhos, means = thermo_rows(folder, manifest)
    z = configurational_z(means, N, rhos)
    dense = folder / 'transport'
    stress = np.atleast_2d(np.loadtxt(dense / 'stress_offdiag.dat'))
    current = np.atleast_2d(np.loadtxt(dense / 'charge_current.dat'))
    com = np.atleast_2d(np.loadtxt(dense / 'global_vcom.dat'))
    assert stress.shape[1] == 6 and current.shape[1] == 4 and com.shape[1] == 4
    assert np.array_equal(stress[:, 0], current[:, 0])
    volume = N / RHO
    sigma = conductivity(current[:, 1:4], 5 * DT, volume)
    np.savetxt(dense / 'conductivity_running.csv',
               np.column_stack([sigma['time'], sigma['correlation'], sigma['running']]),
               delimiter=',', header='lag_tau,current_dot_autocorrelation,sigma_running')
    if smoke:
        # Short smoke trajectories cannot estimate diffusion or a 5-10 tau plateau.
        with __import__('netCDF4').Dataset(dense / 'positions.nc') as dataset:
            position_frames = len(dataset.variables['time'])
        with __import__('netCDF4').Dataset(dense / 'velocities.nc') as dataset:
            velocity_frames = len(dataset.variables['time'])
        assert position_frames >= 2 and velocity_frames >= 2
        result = {'name': manifest['composition']['name'], 'replica': 1,
                  'smoke': True, 'parsed': ['thermo', 'positions', 'velocities',
                  'stress', 'charge_current', 'global_vcom'],
                  'stress_rows': len(stress), 'current_rows': len(current),
                  'position_frames': position_frames, 'velocity_frames': velocity_frames,
                  'sigma_short_cutoffs': sigma['cutoffs']}
        atomic_json(folder / 'smoke_analysis.json', result)
        print('Smoke parsed', folder, flush=True)
        return result
    sex = entropy(rhos, z, means[:, 4], manifest['composition']['q2'])
    assert len(rhos) == 19 and np.isfinite(sex).all()
    time, position, typ, box = read_positions(dense)
    labels = [s['label'] for s in species]
    lag, curve = msd_curve(time, position, typ, len(species))
    fits = fit_msd(lag, curve, WINDOWS)
    assert fits['20_100'] is not None
    blocks = block_diffusion(time, position, typ, labels)
    fluctuations = concentration_fluctuations(position, typ, box, species)
    del position
    diffusion = dict(zip(['net', *labels], fits['20_100']))
    diffusion['probe'] = (diffusion['probe_plus2'] + diffusion['probe_minus2']) / 2
    eta = viscosity(stress[:, 1:6], 5 * DT, volume)
    assert eta['plateau'] is not None
    sigma_ne = nernst_einstein(species, diffusion)
    assert sigma_ne > 0
    sigma_out = {key: value for key, value in sigma.items()
                 if key not in {'time', 'running', 'correlation'}}
    sigma_out['sigma_NE'] = sigma_ne
    sigma_out['inverse_Haven_ratio'] = sigma['plateau'] / sigma_ne
    result = {'name': manifest['composition']['name'], 'code': manifest['composition']['code'],
              'replica': manifest['replica'], 'seed': manifest['seed'], 'N': N,
              's_ex': float(sex[-1]),
              's_ex_pchip': float(entropy(rhos, z, means[:, 4],
                                          manifest['composition']['q2'], True)[-1]),
              'Z': float(z[-1]), 'u': float(means[-1, 4] / T),
              'transport': {'D_net': diffusion['net'], 'D_probe': diffusion['probe'],
                            'D_species': {label: diffusion[label] for label in labels},
                            'fitting_windows': fits, 'labels': ['net', *labels],
                            'block_diffusion': blocks, 'Green_Kubo': vacf(dense, species),
                            'eta': eta, 'D_YH': {key: yeh_hummer(value, eta['plateau'], N)
                                                for key, value in diffusion.items()},
                            'sigma': sigma_out,
                            'concentration_fluctuations': fluctuations,
                            'n_position_frames': len(time)}}
    np.savetxt(folder / 'entropy_density.csv',
               np.column_stack([rhos, z, means[:, 4] / T, sex]), delimiter=',',
               header='rho,Z,u_ex_over_kBT,s_ex_over_kB')
    np.savetxt(dense / 'msd_multi_origin.csv', np.column_stack([lag, curve]),
               delimiter=',', header='lag_tau,net,' + ','.join(labels))
    atomic_json(folder / 'replica_result.json', result)
    print('Analyzed', folder, flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replica', type=Path, required=True)
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    analyze(args.replica, args.smoke)
