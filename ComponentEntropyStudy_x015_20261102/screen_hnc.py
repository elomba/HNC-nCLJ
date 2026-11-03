#!/usr/bin/env python3
"""Reproducible HNC campaign at elevated divalent-probe concentration."""

import argparse
import json
import math
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / "ComponentEntropyStudy_20261030"
BINARY = ROOT.parent / "bin" / "hncr_nclj"
X_PROBE = 0.075
DENSITIES = [0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15,
             0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]
BETAS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def species(label, q, x, eps=1.0):
    return dict(label=label, q=q, x=x, eps=eps, mass=1.0, sigma=1.0)


def candidate(name, entries, description, family):
    entries = [e for e in entries if e["x"] > 1e-12]
    neutral = [e for e in entries if e["q"] == 0]
    charged = [e for e in entries if e["q"] != 0 and not e["label"].startswith("probe_")]
    probes = [e for e in entries if e["label"].startswith("probe_")]
    entries = neutral + charged + probes
    x = np.array([e["x"] for e in entries])
    q = np.array([e["q"] for e in entries])
    assert len(entries) <= 10 and np.all(x > 0)
    assert abs(x.sum() - 1) < 1e-12 and abs(x @ q) < 1e-12
    assert len(probes) == 2 and all(abs(e["x"] - X_PROBE) < 1e-12 for e in probes)
    return dict(name=name, species=entries, description=description, family=family,
                n_components=len(entries), I=float(0.5 * x @ (q * q)),
                ion_fraction=float(x[q != 0].sum()), mixing_entropy=float(-x @ np.log(x)),
                q2=float(x @ (q * q)), q4=float(x @ (q ** 4)),
                mean_sqrt_epsilon=float(sum(e["x"] * math.sqrt(e["eps"]) for e in entries)))


def candidates():
    base = 1 - 2 * X_PROBE
    probes = [species("probe_plus2", 2, X_PROBE), species("probe_minus2", -2, X_PROBE)]

    def add(name, background, description, family):
        return candidate(name, background + probes, description, family)

    cs = [add("baseline_3C", [species("solvent", 0, base)], "Original 3C probe electrolyte", "baseline"),
          add("labels_4C", [species("solvent", 0, base - .30), species("solvent_label", 0, .30)],
              "Mechanically identical neutral-label split", "labels")]
    for f in [.05, .10, .20]:
        cs.append(add(f"mono_f{round(f * 100):03d}",
                      [species("solvent", 0, base - f), species("background_plus1", 1, f / 2),
                       species("background_minus1", -1, f / 2)], "Symmetric monovalent salt", "monovalent"))
    for f in [.04, .10, .20]:
        for t in [0., .25, .5, 1.]:
            cs.append(add(f"asym13_f{round(f * 100):03d}_t{round(t * 100):03d}",
                          [species("solvent", 0, base - f), species("background_plus1", 1, .75 * f * (1 - t)),
                           species("background_minus3", -3, .25 * f * (1 - t)),
                           species("background_minus1", -1, .75 * f * t),
                           species("background_plus3", 3, .25 * f * t)],
                          "Conjugate 1:3 salts", "conjugate_salts"))
    for f in [.10, .30]:
        for eps in [.25, .50, .75]:
            cs.append(add(f"weak_f{round(f * 100):03d}_e{round(eps * 100):03d}",
                          [species("solvent", 0, base - f), species("weak_solvent", 0, f, eps)],
                          "Weak neutral solvent", "weak_solvent"))
    for old in json.loads((OLD / "additional_candidates.json").read_text()):
        entries = [e.copy() for e in old["species"]]
        for e in entries:
            if e["label"] == "solvent":
                e["x"] -= 2 * (X_PROBE - .03)
            elif e["label"].startswith("probe_"):
                e["x"] = X_PROBE
        cs.append(candidate(old["name"], entries, old["description"], old["family"]))
    assert len(cs) == 35 and len({c["name"] for c in cs}) == 35
    by_name = {c["name"]: c for c in cs}
    for group in [("weak_f030_e025", "dispersion_split_d000", "dispersion_split_d025",
                   "dispersion_split_d040"), ("mono_weak_uniform6C", "mono_weak_split7C")]:
        amplitudes = [by_name[name]["mean_sqrt_epsilon"] for name in group]
        assert max(amplitudes) - min(amplitudes) < 1e-12
    return cs


def write_input(path, c, rho, betas, fine=False, restart=False):
    entries = c["species"]
    n = len(entries)
    nr, dr, num = (4096, .005, 100) if fine else (2048, .01, 40)
    lines = [c["name"], f"{nr} {dr} 1.0d-7 0.5 300.0 0.05", f"{num} 1.0d-7", str(n)]
    lines += [" ".join(["1.0"] * n) for _ in entries]
    lines += [" ".join(f"{math.sqrt(a['eps'] * b['eps']):.12f}" for b in entries) for a in entries]
    lines += [" ".join(str(e["q"]) for e in entries), f"1 {len(betas)}", str(rho),
              " ".join(str(b) for b in betas), ".f." if restart else ".t.",
              " ".join(f"{e['x']:.12f}" for e in entries[:-2])]
    lines += [" ".join(["2.5"] * n) for _ in entries] + ["1.0"]
    path.write_text("\n".join(lines) + "\n")


def stability(path, n):
    arr = np.loadtxt(path, ndmin=2)
    assert arr.shape[1] == 1 + n * (n + 1) // 2 + 2
    rows = []
    for row in arr[:20]:
        matrix = np.zeros((n, n))
        ix = 1
        for i in range(n):
            for j in range(i + 1):
                matrix[i, j] = matrix[j, i] = row[ix]
                ix += 1
        eig = np.linalg.eigvalsh(matrix)
        rows.append(dict(k=float(row[0]), min_eigenvalue=float(eig[0]),
                         max_eigenvalue=float(eig[-1]), S_nn=float(row[-2]), S_cc=float(row[-1])))
    return dict(min_eigenvalue_low_k=rows[0]["min_eigenvalue"],
                max_eigenvalue_low_k=rows[0]["max_eigenvalue"],
                max_S_cc_low_k=max(r["S_cc"] for r in rows), low_k=rows)


def solve(c, folder, rho=.7, betas=None, fine=False, restart_from=None):
    folder.mkdir(parents=True, exist_ok=True)
    output = folder / "result.json"
    betas = betas or [.1, .2, .3, .4, .5]
    if output.exists():
        cached = json.loads(output.read_text())
        if (abs(cached["beta"] - betas[-1]) < 1e-10 and abs(cached["rho"] - rho) < 1e-10
                and abs(cached.get("x_probe", .075) - X_PROBE) < 1e-10):
            if cached["status"] == "failed" and not cached.get("error"):
                log_path = folder / "stdout.log"
                detail = log_path.read_text().splitlines()[-12:] if log_path.exists() else []
                cached["error"] = "Solver produced no converged thermodynamic row. " + " ".join(detail)
                output.write_text(json.dumps(cached, indent=2) + "\n")
            return cached
    restart = restart_from is not None and restart_from.exists()
    if restart:
        shutil.copy2(restart_from, folder / "solin.hncr")
    write_input(folder / "input.hncr", c, rho, betas, fine, restart)
    start = time.monotonic()
    try:
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        with (folder / "stdout.log").open("w") as log:
            subprocess.run([str(BINARY), "input.hncr"], cwd=folder, stdout=log,
                           stderr=subprocess.STDOUT, env=env, timeout=600, check=True)
        output_text = (folder / "stdout.log").read_text()
        if "Aborting" in output_text:
            raise RuntimeError("Solver aborted: " + " ".join(output_text.splitlines()[-12:]))
        arr = np.loadtxt(folder / "Thermohncr.dat", comments="#", ndmin=2)
        assert len(arr) == len(betas) and abs(arr[-1, 2] - betas[-1]) < 1e-8
        assert np.isfinite(arr).all()
        vals = arr[-1]
        diagnostic = stability(folder / "skhnc.dat", c["n_components"])
        result = dict(name=c["name"], status="converged", rho=rho, beta=betas[-1], x_probe=X_PROBE,
                      s_ex=float(vals[7]), s2=float(vals[9]), Z=float(vals[3]),
                      chi=float(vals[11]), a_ex=float(vals[6]),
                      u_ex_over_kBT=float(vals[4]), **diagnostic)
    except Exception as exc:
        result = dict(name=c["name"], status="failed", rho=rho, beta=betas[-1],
                      x_probe=X_PROBE, error=str(exc) or "Solver failed without diagnostic")
    result["elapsed_seconds"] = time.monotonic() - start
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(c["name"], folder.parent.name, rho, betas[-1], result["status"], flush=True)
    return result


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    global X_PROBE
    parser = argparse.ArgumentParser()
    parser.add_argument("--x-probe", type=float, default=X_PROBE)
    args = parser.parse_args()
    X_PROBE = args.x_probe
    assert .0 < X_PROBE < .125, "Probe fraction must keep all design solvent fractions positive"
    if not BINARY.exists():
        subprocess.run(["make"], cwd=ROOT.parent, check=True)
    cs = candidates()
    dump(ROOT / "candidates.json", cs)
    by_name = {c["name"]: c for c in cs}
    with ThreadPoolExecutor(max_workers=4) as pool:
        screen = list(pool.map(lambda c: solve(c, ROOT / "hnc" / "screen" / c["name"]), cs))
    dump(ROOT / "hnc_screen_summary.json", screen)
    selected = [c["name"] for c in json.loads((OLD / "selected.json").read_text())]
    best = max(r["s_ex"] for r in screen if r["name"] in selected and r["status"] == "converged")
    extras = [r["name"] for r in screen if r["name"] not in selected and
              r["status"] == "converged" and r["s_ex"] > best]
    fine_names = selected + extras
    def refine(name):
        return solve(by_name[name], ROOT / "hnc" / "fine" / name, fine=True,
                     betas=[.5], restart_from=ROOT / "hnc" / "screen" / name / "lsol.hncr")

    with ThreadPoolExecutor(max_workers=4) as pool:
        fine = list(pool.map(refine, fine_names))
    dump(ROOT / "hnc_fine_summary.json", fine)
    dump(ROOT / "hnc_stability_diagnostics.json",
         [{"name": r["name"], "status": r["status"], "chi": r.get("chi"),
           "min_eigenvalue_low_k": r.get("min_eigenvalue_low_k"),
           "max_eigenvalue_low_k": r.get("max_eigenvalue_low_k"),
           "max_S_cc_low_k": r.get("max_S_cc_low_k")} for r in fine])
    def density_scan(item):
        name, rho = item
        return solve(by_name[name], ROOT / "hnc" / "ti" / name / f"rho_{rho:.4f}", rho=rho)

    with ThreadPoolExecutor(max_workers=4) as pool:
        ti_rows = list(pool.map(density_scan, [(name, rho) for name in selected for rho in DENSITIES]))
    ti = {name: ti_rows[i * len(DENSITIES):(i + 1) * len(DENSITIES)]
          for i, name in enumerate(selected)}

    def temperature_scan(name):
        c = by_name[name]
        rows = []
        previous = None
        for beta in BETAS:
            folder = ROOT / "hnc" / "temperature" / name / f"beta_{beta:.1f}"
            result = solve(c, folder, betas=[beta] if previous else [.1, .2, .3, .4, .5],
                           restart_from=previous)
            rows.append(result)
            if result["status"] == "converged":
                previous = folder / "lsol.hncr"
        return rows

    with ThreadPoolExecutor(max_workers=4) as pool:
        temp_rows = list(pool.map(temperature_scan, selected))
    margins = dict(zip(selected, temp_rows, strict=True))
    dump(ROOT / "hnc_ti_summary.json", ti)
    dump(ROOT / "hnc_temperature_summary.json", margins)
    dump(ROOT / "campaign_manifest.json", dict(x_probe=X_PROBE, screen_attempted=len(screen),
         screen_converged=sum(r["status"] == "converged" for r in screen),
         fine_attempted=len(fine), fine_converged=sum(r["status"] == "converged" for r in fine),
         selected=selected, additional_above_selected=extras))


if __name__ == "__main__":
    main()
