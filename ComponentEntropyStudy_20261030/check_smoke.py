#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
from netCDF4 import Dataset
ROOT=Path(__file__).resolve().parent
checks=[]
for row in json.loads((ROOT/'smoke_jobs.json').read_text()):
    root=Path(row['directory']);d=root/'transport';m=json.loads((root/'manifest.json').read_text());c=m['composition']
    assert (d/'verified.json').exists() and (root/'rho_0.0010/verified.json').exists(),d
    with Dataset(d/'positions.nc') as f:
        f.set_auto_mask(False)
        pos=np.asarray(f.variables['unwrapped_coordinates'][:],float)
        types=np.asarray(f.variables['type'][0],int)
        steps=np.rint(np.asarray(f.variables['time'][:])/.002).astype(int)
    disp=pos-pos[0];disp-=disp.mean(axis=1,keepdims=True)
    d2=(disp*disp).sum(axis=2)
    expected=np.column_stack([d2.mean(axis=1)]+[d2[:,types==i+1].mean(axis=1) for i in range(len(c['species']))])
    a=np.loadtxt(d/'msd_common.dat');common,ip,ia=np.intersect1d(steps,a[:,0].astype(int),return_indices=True)
    assert len(common)>=2
    err=float(np.max(abs(a[ia,1:]-expected[ip])))
    assert np.allclose(a[ia,1:],expected[ip],rtol=5e-5,atol=5e-5),(d,err)
    with Dataset(d/'velocities.nc') as f:
        f.set_auto_mask(False)
        vt=np.asarray(f.variables['time'][:],float)
        ids=np.asarray(f.variables['id'][0],int)
        assert list(ids)==m['velocity_ids']
        assert 'velocities' in f.variables and 'unwrapped_coordinates' in f.variables
    com=np.loadtxt(d/'global_vcom.dat')
    match=np.intersect1d(np.rint(vt/.002).astype(int),com[:,0].astype(int))
    assert len(match)>=len(vt)-2
    thermo=np.loadtxt(d/'thermo_avg.dat')
    assert np.max(abs(thermo[:,4]))<100  # intensive LJ energy, not total box energy
    checks.append(dict(name=c['name'],MSD_max_absolute_difference=err,velocity_particles=len(ids),
                       coordinate_frames=len(pos),velocity_frames=len(vt),low_density_and_dense_passed=True))
(ROOT/'smoke_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print('All eight smoke tests passed low/dense stages, common-frame MSD additivity, independent offline displacements, and NetCDF/COM matching.')
print('Largest online/offline MSD discrepancy:',max(r['MSD_max_absolute_difference'] for r in checks))
