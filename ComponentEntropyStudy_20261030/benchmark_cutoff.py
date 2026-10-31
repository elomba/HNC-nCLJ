#!/usr/bin/env python3
from pathlib import Path
import json,shutil,subprocess
from prepare_md import write_input
ROOT=Path(__file__).resolve().parent
source=ROOT/'smoke/mono_weak_split7C/replica_1'
m=json.loads((source/'manifest.json').read_text())
root=ROOT/'cutoff_benchmark';root.mkdir(exist_ok=False)
folders=[]
for rc in [5,8]:
    d=root/f'coulomb_{rc}'/'transport';d.mkdir(parents=True)
    for name in ['data.mixture','velocity_group.in']:shutil.copy2(source/'transport'/name,d/name)
    (d.parent/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    write_input(d/'in.run',m['composition'],.7,m['seed'],transport=True,smoke=True)
    text=(d/'in.run').read_text().replace('lj/cut/coul/long 2.5 2.5',f'lj/cut/coul/long 2.5 {rc}')
    (d/'in.run').write_text(text)
    folders.append(d)
base=(source/'job.sbatch').read_text().split('cd ',1)[0]
base=base.replace('ce_mono_weak_split7C_r1_smoke','ce_cutoff_benchmark')
base=base.replace(str(source),str(root))
for d in folders:
    base+=f'''cd {d}
mpirun --bind-to none --oversubscribe -np 2 "$HOME/bin/run_lammps_latest" -sf gpu -pk gpu 1 neigh no -in in.run > lammps_run.log 2>&1
python3 {ROOT/'verify_run.py'} {d}
'''
(root/'job.sbatch').write_text(base)
subprocess.run(['sbatch','--test-only',str(root/'job.sbatch')],check=True)
job=subprocess.check_output(['sbatch','--parsable',str(root/'job.sbatch')],text=True).strip()
(root/'submission.json').write_text(json.dumps(dict(job_id=job,purpose='Optimize only the Ewald real-space split at unchanged LJ cutoff and 1e-5 force tolerance'))+'\n')
print('Submitted cutoff benchmark',job,flush=True)
