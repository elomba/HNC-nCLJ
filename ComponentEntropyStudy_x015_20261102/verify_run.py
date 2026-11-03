#!/usr/bin/env python3
"""Verify one EOS or dense LAMMPS stage before advancing an allocation."""
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def verify(folder):
    log = (folder / 'lammps_run.log').read_text(errors='replace')
    assert 'Total wall time:' in log and 'ERROR' not in log, folder
    assert (folder / 'final.restart').is_file() and (folder / 'data.final').is_file(), folder
    assert (folder / 'thermo_avg.dat').is_file(), folder
    result = {'normal_completion': True}
    if folder.name == 'transport':
        manifest = json.loads((folder.parent / 'manifest.json').read_text())
        fractions = np.array([item['x'] for item in manifest['composition']['species']])
        assert abs(sum(item['q'] * item['x'] for item in manifest['composition']['species'])) < 1e-12
        msd = np.atleast_2d(np.loadtxt(folder / 'msd_common.dat'))
        assert np.allclose(msd[:, 1], msd[:, 2:] @ fractions, rtol=5e-5, atol=2e-5)
        stress = np.atleast_2d(np.loadtxt(folder / 'stress_offdiag.dat'))
        current = np.atleast_2d(np.loadtxt(folder / 'charge_current.dat'))
        com = np.atleast_2d(np.loadtxt(folder / 'global_vcom.dat'))
        assert stress.shape[1] == 6 and current.shape[1] == 4 and com.shape[1] == 4
        assert len(stress) == len(current) == len(com)
        assert np.array_equal(stress[:, 0], current[:, 0])
        with Dataset(folder / 'positions.nc') as dataset:
            result['position_frames'] = len(dataset.variables['time'])
        with Dataset(folder / 'velocities.nc') as dataset:
            result['velocity_frames'] = len(dataset.variables['time'])
        expected = manifest['dense_production_steps'] // 500 + 1
        assert abs(result['position_frames'] - expected) <= 1
        assert result['velocity_frames'] >= 2
        result['stress_rows'] = len(stress)
        result['current_rows'] = len(current)
    (folder / 'verified.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Verified', folder, flush=True)


if __name__ == '__main__':
    verify(Path(sys.argv[1]))
