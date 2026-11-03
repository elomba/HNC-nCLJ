#!/usr/bin/env python3
"""Prepare replicated fixed-Hamiltonian EOS and common-frame transport jobs."""
from pathlib import Path
import argparse,hashlib,json,math,os,shlex,subprocess,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
N=4000
DENSITIES=[.001,.0025,.005,.01,.02,.05,.10,.15,.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70]
SEEDS=[41057,52691,67433,83819]
STATES={'baseline_3C':'B3','mono_f010':'B5','weak_f030_e025':'U4',
        'dispersion_split_d000':'L5','dispersion_split_d025':'M5',
        'dispersion_split_d040':'H5','mono_weak_uniform6C':'U6',
        'mono_weak_split7C':'H7'}


def select():
    candidates={c['name']:c for c in json.loads((ROOT/'candidates.json').read_text())}
    selected=[]
    for name,code in STATES.items():
        c=dict(candidates[name]);c['code']=code
        counts=[round(N*s['x']) for s in c['species']]
        assert sum(counts)==N
        assert all(abs(count/N-s['x'])<1e-12 for count,s in zip(counts,c['species']))
        assert sum(count*s['q'] for count,s in zip(counts,c['species']))==0
        assert counts[-2:]==[300,300]
        selected.append(c)
    (ROOT/'selected.json').write_text(json.dumps(selected,indent=2)+'\n')
    return selected


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def write_data(path,c,rho,seed):
    entries=c['species'];counts=[round(N*s['x']) for s in entries]
    assert sum(counts)==N and all(abs(n/N-s['x'])<1e-12 for n,s in zip(counts,entries))
    assert abs(sum(n*s['q'] for n,s in zip(counts,entries)))<1e-10
    rng=np.random.default_rng(seed)
    perm=rng.permutation(N)
    # Reserve identical probe IDs across all compositions of a given replica.
    types=np.zeros(N,dtype=int);cursor=0
    order=[len(entries)-2,len(entries)-1]+list(range(len(entries)-2))
    for index in order:
        ids=perm[cursor:cursor+counts[index]];cursor+=counts[index];types[ids]=index+1
    L=(N/rho)**(1/3);side=math.ceil(N**(1/3));spacing=L/side
    pos=(np.indices((side,side,side)).reshape(3,-1).T[:N]+.5)*spacing
    pos+=(rng.random(pos.shape)-.5)*.06*spacing
    velocity_ids=[]
    for i,sp in enumerate(entries):
        ids=np.flatnonzero(types==i+1)+1
        velocity_ids.extend(ids if sp['q']!=0 else ids[:min(200,len(ids))])
    lines=[f'# {c["name"]}: N={N}, rho={rho}, replica seed={seed}',f'{N} atoms',f'{len(entries)} atom types','',
           f'0.0 {L:.14g} xlo xhi',f'0.0 {L:.14g} ylo yhi',f'0.0 {L:.14g} zlo zhi','','Masses','']
    lines += [f'{i+1} 1.0 # {sp["label"]}' for i,sp in enumerate(entries)]+['','Atoms # charge','']
    lines += [f'{i+1} {types[i]} {entries[types[i]-1]["q"]} '+ ' '.join(f'{v:.14g}' for v in pos[i]) for i in range(N)]
    path.write_text('\n'.join(lines)+'\n')
    return sorted(int(i) for i in velocity_ids)


def write_input(path,c,rho,seed,transport=False,smoke=False):
    entries=c['species'];n=len(entries)
    # Fixed fractional-box Ewald real-space split: 5 sigma at rho=.7.
    # Only numerical splitting changes; every LJ pair cutoff stays 2.5.
    rc=5*(.7/rho)**(1/3)
    assert rc<.5*(N/rho)**(1/3)
    eq=100000 if transport else 50000
    prod=1500000 if transport else 100000
    if smoke:eq,prod=20,2000 if transport else 1000
    lines=[f'# {c["name"]}: fixed LJ cutoff 2.5 at every density',
           f'variable n_eq index {eq}',f'variable n_prod index {prod}',
           'units lj','atom_style charge','boundary p p p','read_data data.mixture',
           f'pair_style lj/cut/coul/long 2.5 {rc:.14g}']
    for i,a in enumerate(entries):
        for jj in range(i,n):
            b=entries[jj];lines.append(f'pair_coeff {i+1} {jj+1} {math.sqrt(a["eps"]*b["eps"]):.14g} 1.0 2.5')
    lines+=['pair_modify shift yes tail no','kspace_style pppm 1.0e-5','timestep 0.002',
            'neighbor 0.3 bin','neigh_modify delay 0 every 1 check yes',
            'thermo 1000','thermo_style custom step temp press pe ke etotal density',
            'min_style cg','minimize 1e-4 1e-6 2000 10000',
            f'velocity all create 2.0 {seed} dist gaussian','velocity all zero linear',
            'fix integrator all nvt temp 2.0 2.0 0.1','run ${n_eq}','reset_timestep 0',
            'variable temp0 equal temp','variable press0 equal press','variable pe0 equal pe',
            f'variable z0 equal press/({rho:.14g}*2.0)',
            'fix thermo_output all ave/time 10 50 500 v_temp0 v_press0 v_z0 v_pe0 file thermo_avg.dat',
            f'restart {250000 if transport else 100000} checkpoint.*.restart']
    if transport:
        # Common global COM correction; no group-specific drift subtraction.
        for i in range(n):lines.append(f'group species{i+1} type {i+1}')
        lines+=['compute gcom all com','compute displacement all displace/atom',
                'thermo_style custom step temp press pe ke etotal density c_gcom[1] c_gcom[2] c_gcom[3]',
                'run 0 post no',
                'variable com0x equal $(c_gcom[1]:%.16g)',
                'variable com0y equal $(c_gcom[2]:%.16g)',
                'variable com0z equal $(c_gcom[3]:%.16g)',
                'variable gdx equal c_gcom[1]-v_com0x',
                'variable gdy equal c_gcom[2]-v_com0y',
                'variable gdz equal c_gcom[3]-v_com0z',
                'variable dx atom c_displacement[1]-v_gdx',
                'variable dy atom c_displacement[2]-v_gdy',
                'variable dz atom c_displacement[3]-v_gdz',
                'variable r2 atom v_dx*v_dx+v_dy*v_dy+v_dz*v_dz',
                'compute msd_all all reduce ave v_r2']
        for i in range(n):lines.append(f'compute msd_{i+1} species{i+1} reduce ave v_r2')
        lines.append('fix msd_output all ave/time 100 1 100 c_msd_all '+' '.join(f'c_msd_{i+1}' for i in range(n))+' file msd_common.dat')
        lines+=['dump position_output all netcdf 500 positions.nc id type xu yu zu',
                'include velocity_group.in',
                'variable vcomx equal vcm(all,x)','variable vcomy equal vcm(all,y)','variable vcomz equal vcm(all,z)']
        lines+=['variable sxy equal pxy','variable sxz equal pxz','variable syz equal pyz',
                'variable sd1 equal (pxx-pyy)/2','variable sd2 equal (pyy-pzz)/2',
                'fix stress_output all ave/time 5 1 5 v_sxy v_sxz v_syz v_sd1 v_sd2 file stress_offdiag.dat',
                'variable jx atom q*vx','variable jy atom q*vy','variable jz atom q*vz',
                'compute charge_current all reduce sum v_jx v_jy v_jz',
                'fix charge_output all ave/time 5 1 5 c_charge_current[1] c_charge_current[2] c_charge_current[3] file charge_current.dat',
                'fix vcom_output all ave/time 5 1 5 v_vcomx v_vcomy v_vcomz file global_vcom.dat']
        if not smoke:lines.append('run 1450000')
        else:lines.append('run 1500')
        lines+=['dump velocity_output velocity_sample netcdf 5 velocities.nc id type xu yu zu vx vy vz',
                'run 500'] if smoke else ['dump velocity_output velocity_sample netcdf 5 velocities.nc id type xu yu zu vx vy vz',
                'run 50000']
    else:lines.append('run ${n_prod}')
    lines+=['write_data data.final pair ij','write_restart final.restart']
    path.write_text('\n'.join(lines)+'\n')


def prepare(selected,smoke=False):
    root=ROOT/('smoke' if smoke else 'md')
    root.mkdir(exist_ok=True)
    jobs=[]
    for c in selected:
        for replica,seed in enumerate(SEEDS if not smoke else SEEDS[:1],1):
            folder=root/c['name']/f'replica_{replica}'
            folder.mkdir(parents=True,exist_ok=False)
            runs=[]
            densities=([.001] if smoke else DENSITIES[:-1])
            for rho in densities:
                d=folder/f'rho_{rho:.4f}';d.mkdir()
                ids=write_data(d/'data.mixture',c,rho,seed)
                write_input(d/'in.run',c,rho,seed,smoke=smoke)
                runs.append(str(d))
            transport=folder/'transport';transport.mkdir()
            ids=write_data(transport/'data.mixture',c,.7,seed)
            (transport/'velocity_group.in').write_text('group velocity_sample id '+' '.join(str(i) for i in ids)+'\n')
            write_input(transport/'in.run',c,.7,seed,transport=True,smoke=smoke)
            runs.append(str(transport))
            manifest=dict(composition=c,smoke=smoke,N=N,T=2.,rho_transport=.7,seed=seed,replica=replica,
                          densities=([.001,.7] if smoke else DENSITIES),velocity_ids=ids,runs=runs,
                          dense_equilibration_steps=20 if smoke else 100000,
                          dense_production_steps=2000 if smoke else 1500000,
                          dense_velocity_steps=500 if smoke else 50000,
                          integrator='NVT, Tdamp=0.1, dt=0.002',
                          frame='one global center of mass for every species',
                          input_hashes={str(Path(r)/'in.run'):digest(Path(r)/'in.run') for r in runs},
                          data_hashes={str(Path(r)/'data.mixture'):digest(Path(r)/'data.mixture') for r in runs})
            (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            name=f'x15_{c["code"]}_r{replica}'+('_smoke' if smoke else '')
            content=f'''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --partition=gpu_v4_short
#SBATCH --nodelist=ladon28,ladon29
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:nvidia_rtx_pro_4500_blackwell:1
#SBATCH --mem=8G
#SBATCH --time={'00:30:00' if smoke else '12:00:00'}
#SBATCH --output={folder}/slurm_%j.out
#SBATCH --error={folder}/slurm_%j.err
set -eo pipefail
source "$HOME/bin/setup_trj"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_coll=^hcoll
'''
            for r in runs:
                content+=f'''cd {shlex.quote(r)}
echo "Starting {Path(r).name} on $HOSTNAME at $(date)"
if test -e final.restart || compgen -G 'checkpoint.*.restart' > /dev/null; then
    echo 'Existing output/checkpoint detected; dense-stage resume is unsupported' >&2
    exit 2
fi
mpirun --bind-to none --oversubscribe -np 2 "$HOME/bin/run_lammps_latest" -sf gpu -pk gpu 1 neigh no -in in.run > lammps_run.log 2>&1
python3 {shlex.quote(str(ROOT/'verify_run.py'))} {shlex.quote(r)}
'''
            content+=f'python3 {shlex.quote(str(ROOT/"analyze_md.py"))} --replica {shlex.quote(str(folder))}'+(' --smoke\n' if smoke else '\n')
            content+=f'echo "Completed {name} at $(date)"\n'
            (folder/'job.sbatch').write_text(content)
            jobs.append(dict(name=name,directory=str(folder),composition=c['name'],replica=replica))
    (ROOT/('smoke_jobs.json' if smoke else 'jobs.json')).write_text(json.dumps(jobs,indent=2)+'\n')
    print('Prepared',len(jobs),'GPU allocation scripts;',len(jobs)*(2 if smoke else len(DENSITIES)),'simulation stages')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--smoke',action='store_true');args=p.parse_args()
    selection=select()
    prepare(selection,args.smoke)
