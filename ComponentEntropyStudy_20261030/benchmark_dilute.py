#!/usr/bin/env python3
from pathlib import Path
import json,shutil,subprocess
ROOT=Path(__file__).resolve().parent
source=ROOT/'smoke/mono_weak_split7C/replica_1';d=ROOT/'dilute_cutoff_benchmark/rho_0.0010';d.mkdir(parents=True,exist_ok=False)
shutil.copy2(source/'rho_0.0010/data.mixture',d/'data.mixture')
rc=5*(.7/.001)**(1/3)
text=(source/'rho_0.0010/in.run').read_text().replace('lj/cut/coul/long 2.5 8.0',f'lj/cut/coul/long 2.5 {rc:.14g}')
(d/'in.run').write_text(text)
header=(source/'job.sbatch').read_text().split('cd ',1)[0]
header=header.replace('ce_mono_weak_split7C_r1_smoke','ce_dilute_cutoff_benchmark').replace(str(source),str(d.parent))
header+=f'''cd {d}
mpirun --bind-to none --oversubscribe -np 2 "$HOME/bin/run_lammps_latest" -sf gpu -pk gpu 1 neigh no -in in.run > lammps_run.log 2>&1
python3 {ROOT/'verify_run.py'} {d}
'''
(d.parent/'job.sbatch').write_text(header)
subprocess.run(['sbatch','--test-only',str(d.parent/'job.sbatch')],check=True)
job=subprocess.check_output(['sbatch','--parsable',str(d.parent/'job.sbatch')],text=True).strip()
(d.parent/'submission.json').write_text(json.dumps(dict(job_id=job,coulomb_cutoff=rc,lj_cutoff=2.5,force_tolerance=1e-5))+'\n')
print('Submitted dilute cutoff benchmark',job,flush=True)
