#!/usr/bin/env python3
"""Summarize and plot the high-divalent HNC campaign."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / "ComponentEntropyStudy_20261030"
RESULTS = ROOT / "results"
LABELS = {"baseline_3C": "B3", "mono_f010": "B5", "weak_f030_e025": "U4",
          "dispersion_split_d000": "L5", "dispersion_split_d025": "M5",
          "dispersion_split_d040": "H5", "mono_weak_uniform6C": "U6",
          "mono_weak_split7C": "H7"}
CONTRASTS = [("H5", "U4"), ("M5", "U4"), ("L5", "U4"),
             ("H7", "U6"), ("H5", "B3"), ("H7", "B5"), ("B5", "B3")]


def read(name, old=False):
    return json.loads(((OLD if old else ROOT) / name).read_text())


def indexed(rows):
    return {r["name"]: r for r in rows}


def flags(row):
    if row["status"] != "converged":
        return ["solver failed"]
    found = []
    if row["chi"] <= 0:
        found.append("chi <= 0")
    if row["min_eigenvalue_low_k"] <= 1e-4:
        found.append("low-k eigenvalue <= 1e-4")
    if row["max_S_cc_low_k"] >= 1:
        found.append("S_cc >= 1")
    return found


def main():
    RESULTS.mkdir(exist_ok=True)
    manifest = read("campaign_manifest.json")
    candidates = read("candidates.json")
    old_screen = indexed(read("hnc_screen_summary.json", old=True))
    new_screen = indexed(read("hnc_screen_summary.json"))
    old_fine = indexed(read("hnc_fine_summary.json", old=True))
    new_fine = indexed(read("hnc_fine_summary.json"))
    ti = read("hnc_ti_summary.json")
    temp = read("hnc_temperature_summary.json")
    names = [c["name"] for c in candidates]
    selected = manifest["selected"]
    comparison = {}
    for name in names:
        old = old_fine.get(name, old_screen[name])
        new = new_fine.get(name, new_screen[name])
        comparison[name] = {"old": {k: old.get(k) for k in ("status", "s_ex", "s2", "Z", "chi")},
                            "new": {k: new.get(k) for k in ("status", "s_ex", "s2", "Z", "chi")},
                            "basis": "fine" if name in old_fine and name in new_fine else "coarse"}
    reverse = {v: k for k, v in LABELS.items()}
    contrasts = {}
    for a, b in CONTRASTS:
        aa, bb = reverse[a], reverse[b]
        contrasts[f"{a}-{b}"] = {
            era: (round(comparison[aa][era]["s_ex"] - comparison[bb][era]["s_ex"], 4)
                  if comparison[aa][era]["s_ex"] is not None and comparison[bb][era]["s_ex"] is not None else None)
            for era in ("old", "new")}
    assert contrasts["L5-U4"]["old"] == contrasts["L5-U4"]["new"] == 0.0
    ti_flags = {name: [{"rho": r["rho"], "flags": flags(r)} for r in rows if flags(r)]
                for name, rows in ti.items()}
    temp_flags = {name: [{"beta": r["beta"], "flags": flags(r)} for r in rows if flags(r)]
                  for name, rows in temp.items()}
    lowest_t = {name: min((1 / r["beta"] for r in rows if r["status"] == "converged"
                           and r["min_eigenvalue_low_k"] > 0), default=None)
                for name, rows in temp.items()}
    summary = dict(manifest=manifest, comparison=comparison, contrasts=contrasts,
                   ti=ti, temperature=temp, ti_flags=ti_flags, temperature_flags=temp_flags,
                   lowest_positive_converged_temperature=lowest_t)
    (RESULTS / "hnc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "figure.dpi": 150})
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), constrained_layout=True)
    xx = np.arange(len(names))
    for offset, era, marker in [(-.16, "old", "o"), (.16, "new", "s")]:
        yy = [comparison[name][era]["s_ex"] for name in names]
        axes[0].scatter(xx + offset, yy, s=25, marker=marker,
                        label=f"x(+2)=x(-2)={.03 if era == 'old' else .075:g}")
    axes[0].set_xticks(xx, [LABELS.get(n, n) for n in names], rotation=75, ha="right")
    axes[0].set_ylabel(r"$s_{\rm ex}/Nk_B$ at $\rho^*=0.7$, $T^*=2$")
    axes[0].legend()
    axes[0].grid(alpha=.25)
    for name in selected:
        rows = [r for r in ti[name] if r["status"] == "converged"]
        axes[1].plot([r["rho"] for r in rows], [r["s_ex"] for r in rows],
                     marker=".", lw=1.3, label=LABELS[name])
    axes[1].set(xlabel=r"Density $\rho^*$", ylabel=r"$s_{\rm ex}/Nk_B$ at $T^*=2$")
    axes[1].legend(ncol=4)
    axes[1].grid(alpha=.25)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"hnc_screening.{ext}")
    plt.close(fig)

    lines = ["# High-divalent HNC report", "", "At T*=2, rho*=0.7; x(+2)=x(-2) changes from 0.03 to 0.075. "
             "The solver's chi column is reduced inverse isothermal compressibility.", "",
             f"Coarse: {manifest['screen_attempted']} attempted, {manifest['screen_converged']} converged. "
             f"Fine: {manifest['fine_attempted']} attempted, {manifest['fine_converged']} converged.", "",
             "Additional candidates exceeding the best selected coarse s_ex: " +
             (", ".join(manifest["additional_above_selected"]) or "none") + ".", "",
             "## All-candidate comparison", "",
             "The eight selected analogues use fine values; other candidates use coarse values.", "",
             "| Candidate | old s_ex | new s_ex | old s2 | new s2 | old Z | new Z | old chi | new chi |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name in names:
        old, new = comparison[name]["old"], comparison[name]["new"]
        def fmt(row, key):
            return f"{row[key]:.4f}" if row[key] is not None else "failed"
        lines.append("| " + name + " | " + " | ".join(
            [fmt(old, "s_ex"), fmt(new, "s_ex"), fmt(old, "s2"), fmt(new, "s2"),
             fmt(old, "Z"), fmt(new, "Z"), fmt(old, "chi"), fmt(new, "chi")]) + " |")
    lines += ["", "## Excess-entropy contrasts", "", "| Contrast | old | new |", "|---|---:|---:|"]
    for label, item in contrasts.items():
        lines.append(f"| {label} | {item['old']:+.4f} | {item['new']:+.4f} |")
    lines += ["", "## TI-path diagnostics", "",
              "Flags: solver failure, chi <= 0, first-k eigenvalue <= 1e-4, or maximum S_cc over first 20 k >= 1. "
              "The small charge-screening eigenvalue at finite k needs physical interpretation; the threshold is a screening alert.", "",
              "| Mixture | rho | status | Z | s_ex | chi | min eigenvalue | S_cc(k1) | flag |",
              "|---|---:|---|---:|---:|---:|---:|---:|---|"]
    for name in selected:
        for r in ti[name]:
            lines.append(f"| {LABELS[name]} | {r['rho']:.4f} | {r['status']} | " +
                         (f"{r['Z']:.4f} | {r['s_ex']:.4f} | {r['chi']:.4f} | "
                          f"{r['min_eigenvalue_low_k']:.6g} | {r['low_k'][0]['S_cc']:.6g} | "
                          if r["status"] == "converged" else "— | — | — | — | — | ") +
                         (", ".join(flags(r)) or "—") + " |")
    lines += ["", "## Temperature margin at rho*=0.7", "",
              "This continuation only probes the HNC closure: failures at lower T* reflect the breakdown of HNC under strong divalent Coulomb coupling, not a phase boundary. The MD campaign stays at T*=2, as in the original study.", "",
              "| Mixture | beta | T* | status | min eigenvalue | S_cc(k1) | chi | flag |",
              "|---|---:|---:|---|---:|---:|---:|---|"]
    for name in selected:
        for r in temp[name]:
            lines.append(f"| {LABELS[name]} | {r['beta']:.1f} | {1/r['beta']:.3f} | {r['status']} | " +
                         (f"{r['min_eigenvalue_low_k']:.6g} | {r['low_k'][0]['S_cc']:.6g} | "
                          f"{r['chi']:.4f} | " if r["status"] == "converged" else "— | — | — | ") +
                         (", ".join(flags(r)) or "—") + " |")
    lines += ["", "## Limitations", "", "HNC convergence and finite-k positivity do not establish phase equilibrium. "
              "The TI path requires MD validation before thermodynamic claims.", ""]
    (RESULTS / "hnc_report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
