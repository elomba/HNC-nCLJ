#!/usr/bin/env python3
"""Generate immutable dense transport inputs for both follow-up campaigns."""
import argparse
import hashlib
import json
import math
import shlex
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent
STATES = {
    "B3": "baseline_3C",
    "U4": "weak_f030_e025",
    "H5": "dispersion_split_d040",
    "U6": "mono_weak_uniform6C",
    "H7": "mono_weak_split7C",
}
CAMPAIGNS = {
    "long_N4000": (4000, [118543, 230917, 347221, 459803, 562129, 673351]),
    "N8000": (8000, [784913, 812477, 925691, 137359]),
}
RHO = 0.7
DT = 0.002


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_data(path, composition, n_atoms, seed):
    species = composition["species"]
    counts = [round(n_atoms * item["x"]) for item in species]
    assert sum(counts) == n_atoms
    assert all(abs(count / n_atoms - item["x"]) < 1e-12 for count, item in zip(counts, species))
    assert sum(count * item["q"] for count, item in zip(counts, species)) == 0
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_atoms)
    types = np.zeros(n_atoms, dtype=int)
    cursor = 0
    for index in [len(species) - 2, len(species) - 1, *range(len(species) - 2)]:
        ids = perm[cursor:cursor + counts[index]]
        cursor += counts[index]
        types[ids] = index + 1
    length = (n_atoms / RHO) ** (1 / 3)
    side = math.ceil(n_atoms ** (1 / 3))
    spacing = length / side
    positions = (np.indices((side, side, side)).reshape(3, -1).T[:n_atoms] + 0.5) * spacing
    positions += (rng.random(positions.shape) - 0.5) * 0.06 * spacing
    ion_ids = sorted(int(i + 1) for i in np.flatnonzero(
        np.array([species[t - 1]["q"] != 0 for t in types])
    ))
    lines = [
        f'# {composition["name"]}: N={n_atoms}, rho={RHO}, seed={seed}',
        f"{n_atoms} atoms", f"{len(species)} atom types", "",
        f"0.0 {length:.14g} xlo xhi", f"0.0 {length:.14g} ylo yhi",
        f"0.0 {length:.14g} zlo zhi", "", "Masses", "",
    ]
    lines += [f'{i + 1} 1.0 # {item["label"]}' for i, item in enumerate(species)]
    lines += ["", "Atoms # charge", ""]
    lines += [
        f'{i + 1} {types[i]} {species[types[i] - 1]["q"]} '
        + " ".join(f"{value:.14g}" for value in positions[i])
        for i in range(n_atoms)
    ]
    path.write_text("\n".join(lines) + "\n")
    return ion_ids


def write_input(path, composition, n_atoms, seed, eq_steps, production_steps, velocity_steps):
    species = composition["species"]
    assert 5 < (n_atoms / RHO) ** (1 / 3) / 2
    assert production_steps > velocity_steps
    lines = [
        f'# {composition["name"]}: dense transport, N={n_atoms}',
        "units lj", "atom_style charge", "boundary p p p", "read_data data.mixture",
        "pair_style lj/cut/coul/long 2.5 5.0",
    ]
    for i, item in enumerate(species):
        for j in range(i, len(species)):
            other = species[j]
            lines.append(
                f'pair_coeff {i + 1} {j + 1} {math.sqrt(item["eps"] * other["eps"]):.14g} 1.0 2.5'
            )
    lines += [
        "pair_modify shift yes tail no", "kspace_style pppm 1.0e-5",
        "timestep 0.002", "neighbor 0.3 bin",
        "neigh_modify delay 0 every 1 check yes", "thermo 1000",
        "thermo_style custom step temp press pe ke etotal density",
        "min_style cg", "minimize 1e-4 1e-6 2000 10000",
        f"velocity all create 2.0 {seed} dist gaussian", "velocity all zero linear",
        "fix integrator all nvt temp 2.0 2.0 0.1",
        f"run {eq_steps}", "reset_timestep 0",
        "variable temp0 equal temp", "variable press0 equal press",
        "variable pe0 equal pe", f"variable z0 equal press/({RHO:.14g}*2.0)",
        "fix thermo_output all ave/time 10 50 500 v_temp0 v_press0 v_z0 v_pe0 file thermo_avg.dat",
        "restart 250000 checkpoint.*.restart",
    ]
    lines += [f"group species{i + 1} type {i + 1}" for i in range(len(species))]
    lines += [
        "compute gcom all com", "compute displacement all displace/atom",
        "thermo_style custom step temp press pe ke etotal density c_gcom[1] c_gcom[2] c_gcom[3]",
        "run 0 post no", "variable com0x equal $(c_gcom[1]:%.16g)",
        "variable com0y equal $(c_gcom[2]:%.16g)",
        "variable com0z equal $(c_gcom[3]:%.16g)",
        "variable gdx equal c_gcom[1]-v_com0x",
        "variable gdy equal c_gcom[2]-v_com0y",
        "variable gdz equal c_gcom[3]-v_com0z",
        "variable dx atom c_displacement[1]-v_gdx",
        "variable dy atom c_displacement[2]-v_gdy",
        "variable dz atom c_displacement[3]-v_gdz",
        "variable r2 atom v_dx*v_dx+v_dy*v_dy+v_dz*v_dz",
        "compute msd_all all reduce ave v_r2",
    ]
    lines += [f"compute msd_{i + 1} species{i + 1} reduce ave v_r2" for i in range(len(species))]
    lines += [
        "fix msd_output all ave/time 500 1 500 c_msd_all "
        + " ".join(f"c_msd_{i + 1}" for i in range(len(species)))
        + " file msd_common.dat",
        "dump position_output all netcdf 500 positions.nc id type xu yu zu",
        "fix gcom_output all ave/time 500 1 500 c_gcom[1] c_gcom[2] c_gcom[3] file global_com.dat",
        "variable sxy equal pxy", "variable sxz equal pxz", "variable syz equal pyz",
        "variable sd1 equal (pxx-pyy)/2", "variable sd2 equal (pyy-pzz)/2",
        "fix stress_output all ave/time 5 1 5 v_sxy v_sxz v_syz v_sd1 v_sd2 file stress_offdiag.dat",
        f"run {production_steps - velocity_steps}",
        "include velocity_group.in",
        "variable vcomx equal vcm(all,x)", "variable vcomy equal vcm(all,y)",
        "variable vcomz equal vcm(all,z)",
        "dump velocity_output velocity_sample netcdf 5 velocities.nc id type xu yu zu vx vy vz",
        "fix vcom_output all ave/time 5 1 5 v_vcomx v_vcomy v_vcomz file global_vcom.dat",
        f"run {velocity_steps}", "write_data data.final pair ij", "write_restart final.restart",
    ]
    path.write_text("\n".join(lines) + "\n")


def prepare(smoke=False):
    selected = {item["name"]: item for item in json.loads((OLD / "selected.json").read_text())}
    root = ROOT / ("smoke" if smoke else "md")
    jobs = []
    for campaign, (n_atoms, seeds) in CAMPAIGNS.items():
        for code, name in STATES.items():
            composition = selected[name]
            for replica, seed in enumerate(seeds[:1] if smoke else seeds, 1):
                folder = root / campaign / name / f"replica_{replica}"
                folder.mkdir(parents=True, exist_ok=False)
                eq_steps = 20 if smoke else 100000
                production_steps = 2000 if smoke else 2500000
                velocity_steps = 500 if smoke else 100000
                ion_ids = write_data(folder / "data.mixture", composition, n_atoms, seed)
                (folder / "velocity_group.in").write_text(
                    "group velocity_sample id " + " ".join(map(str, ion_ids)) + "\n"
                )
                write_input(folder / "in.run", composition, n_atoms, seed, eq_steps,
                            production_steps, velocity_steps)
                manifest = {
                    "composition": composition, "campaign": campaign, "N": n_atoms,
                    "rho": RHO, "T": 2.0, "seed": seed, "replica": replica,
                    "smoke": smoke, "equilibration_steps": eq_steps,
                    "production_steps": production_steps, "velocity_steps": velocity_steps,
                    "input_hashes": {"in.run": digest(folder / "in.run")},
                    "data_hashes": {"data.mixture": digest(folder / "data.mixture")},
                }
                (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
                job_name = f"cf_{campaign}_{code}_r{replica}" + ("_smoke" if smoke else "")
                script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=gpu_v4_short
#SBATCH --nodelist=ladon28,ladon29
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:nvidia_rtx_pro_4500_blackwell:1
#SBATCH --mem=16G
#SBATCH --time={"00:30:00" if smoke else "12:00:00"}
#SBATCH --output={folder}/slurm_%j.out
#SBATCH --error={folder}/slurm_%j.err
set -eo pipefail
source "$HOME/bin/setup_trj"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_coll=^hcoll
cd {shlex.quote(str(folder))}
mpirun --bind-to none --oversubscribe -np 2 "$HOME/bin/run_lammps_latest" -sf gpu -pk gpu 1 neigh no -in in.run > lammps_run.log 2>&1
python3 {shlex.quote(str(ROOT / "verify_followup.py"))} {shlex.quote(str(folder))}
'''
                if not smoke:
                    script += (f'python3 {shlex.quote(str(ROOT / "analyze_followup.py"))} '
                               f'--replica {shlex.quote(str(folder))}\n')
                (folder / "job.sbatch").write_text(script)
                jobs.append({
                    "name": job_name, "campaign": campaign, "directory": str(folder),
                    "composition": name, "replica": replica, "N": n_atoms,
                })
    target = ROOT / ("smoke_jobs.json" if smoke else "jobs.json")
    target.write_text(json.dumps(jobs, indent=2) + "\n")
    print(f"Prepared {len(jobs)} {'smoke' if smoke else 'production'} replicas")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    prepare(parser.parse_args().smoke)
