#!/usr/bin/env python3
"""Synthetic checks for the high-divalent MD estimators."""
import numpy as np
import subprocess
import tempfile
from pathlib import Path

from analyze_md import (charge_current, conductivity, conductivity_integral,
                        entropy, nernst_einstein)
from summarize_md import metric, pooled_inverse_haven, ratio_contrast


def main():
    # Deterministic exponential correlation checks the GK integration formula.
    dt = .001
    lag = np.arange(30001) * dt
    correlation = 6 * np.exp(-lag / 2)
    integral = conductivity_integral(correlation, dt, 10, 2)
    assert abs(integral[-1] - .2) < 1e-6
    # Also exercise the FFT path with a known zero-current series.
    assert abs(conductivity(np.zeros((4000, 3)), .01, 10)['plateau']) < 1e-12
    species = [{'label': 'plus', 'x': .075, 'q': 2},
               {'label': 'minus', 'x': .075, 'q': -2}]
    assert np.isclose(nernst_einstein(species, {'plus': .2, 'minus': .3}), .0525)
    velocities = np.arange(18, dtype=float).reshape(3, 2, 3)
    q = [2, -2]
    assert np.allclose(charge_current(q, velocities),
                       charge_current(q, velocities + np.array([7., -3., 9.])))
    paired = ratio_contrast([2., 3., 4.], [1., 2., 2.])
    assert np.isclose(paired['percent'], 100 * np.expm1(np.mean(np.log([2, 1.5, 2]))))
    synthetic = {'transport': {'sigma': {'cutoffs': {'2': .12, '5': .15,
                                                    '10': .17, '20': .18},
                                         'sigma_NE': .25, 'inverse_Haven_ratio': .72}}}
    for lag, value in [('2', .12), ('5', .15), ('10', .17), ('20', .18)]:
        assert metric(synthetic, f'sigma_{lag}') == value
    assert metric(synthetic, 'sigma_NE') == .25
    assert metric(synthetic, 'inv_Haven') == .72
    pooled = pooled_inverse_haven({'B3': {1: synthetic},
                                   'L5': {1: {'transport': {'sigma': {'inverse_Haven_ratio': 9}}}},
                                   'U4': {1: synthetic}})
    assert pooled['n'] == 2 and pooled['mean'] == .72
    rho = np.array([.001, .01, .1, .7])
    q2 = .6
    c = (4*np.pi*.5*q2)**1.5/(12*np.pi)
    # For Z-1=-C*sqrt(rho)/2 and PE=0, the TI integral is -C*sqrt(rho).
    z = 1 - c*np.sqrt(rho)/2
    assert np.allclose(entropy(rho, z, np.zeros_like(rho), q2), c*np.sqrt(rho))
    job = Path(__file__).resolve().parent / 'smoke/baseline_3C/replica_1/job.sbatch'
    script = job.read_text()
    guard = script[script.index('if test -e final.restart'):]
    guard = guard[:guard.index('\nfi') + 3]
    with tempfile.TemporaryDirectory() as directory:
        assert subprocess.run(['bash', '-c', guard], cwd=directory, check=False,
                              capture_output=True).returncode == 0
        (Path(directory) / 'checkpoint.250000.restart').touch()
        assert subprocess.run(['bash', '-c', guard], cwd=directory, check=False,
                              capture_output=True).returncode == 2
    print('MD analysis synthetic tests passed')


if __name__ == '__main__':
    main()
