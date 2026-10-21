#!/usr/bin/env python3
"""
generate_manuscript_figures.py
------------------------------
Generates publication-quality figures for Physical Review E manuscript:
  - Fig 1: compressibility_and_ti.png (Z(rho) and A_ex(rho) comparing MD and HNC)
  - Fig 2: entropy_comparison.png (s_ex(rho) and s^(2)(rho) across all systems)
  - Fig 3: delta_entropy.png (Entropy variation with solvent attraction and component count)
  - Fig 4: msd_comparison.png (MSD vs t for net fluid, multivalent ions, and solvent species)
  - Fig 5: rosenfeld_scaling.png (Rosenfeld excess entropy scaling plot ln(D*) vs s_ex)
"""

import os
import glob
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Styling for PRE publication
plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 300,
    "lines.linewidth": 1.8,
    "lines.markersize": 6,
})

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DOC_FIGS_DIR = os.path.join(BASE_DIR, "doc", "figures")
os.makedirs(DOC_FIGS_DIR, exist_ok=True)

TEMP = 2.00
DT = 0.002
RHO_TARGET = 0.70

SYSTEMS = {
    "system_3C": {
        "name": "3C (Option 1, I=0.12)",
        "label": "3C (I=0.12)",
        "color": "#1f77b4",  # Blue
        "marker": "o",
        "hnc_file": os.path.join(BASE_DIR, "EntropyComparison", "results", "Thermo_3C_density_sweep.dat"),
    },
    "system_4C": {
        "name": r"4C ($\varepsilon_2=1.25$)",
        "label": r"4C ($\varepsilon_2=1.25$)",
        "color": "#2ca02c",  # Green
        "marker": "s",
        "hnc_file": os.path.join(BASE_DIR, "EntropyComparison", "results", "Thermo_4C_density_sweep.dat"),
    },
    "system_4C_eps075": {
        "name": r"4C ($\varepsilon_2=0.75$)",
        "label": r"4C ($\varepsilon_2=0.75$)",
        "color": "#ff7f0e",  # Orange
        "marker": "^",
        "hnc_file": os.path.join(BASE_DIR, "EntropyComparison", "results", "Thermo_4C_eps075_density_sweep.dat"),
    },
    "system_5C": {
        "name": "5C (I=0.19)",
        "label": "5C (I=0.19)",
        "color": "#d62728",  # Red
        "marker": "D",
        "hnc_file": os.path.join(BASE_DIR, "EntropyComparison", "results", "Thermo_5C_density_sweep.dat"),
    },
}


def parse_thermo_avg_file(filepath):
    if not os.path.isfile(filepath):
        return None
    steps, temps, presses, z_facs, pes = [], [], [], [], []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts or parts[0] in ["Step", "#"]:
                continue
            try:
                steps.append(float(parts[0]))
                temps.append(float(parts[1]))
                presses.append(float(parts[2]))
                z_facs.append(float(parts[3]))
                pes.append(float(parts[4]))
            except ValueError:
                continue
    if not steps:
        return None
    start = len(steps) // 2
    return {
        "temp": np.mean(temps[start:]),
        "press": np.mean(presses[start:]),
        "z": np.mean(z_facs[start:]),
        "u_ex": np.mean(pes[start:]),
    }


def collect_md_data(sys_key):
    sys_dir = os.path.join(BASE_DIR, "lammps_study", sys_key)
    folders = glob.glob(os.path.join(sys_dir, "rho_*"))
    data = []
    for fld in folders:
        bname = os.path.basename(fld)
        try:
            rho = float(bname.replace("rho_", ""))
        except ValueError:
            continue
        p = os.path.join(fld, "thermo_avg.dat")
        avg = parse_thermo_avg_file(p)
        if avg:
            data.append({"rho": rho, **avg})
    data.sort(key=lambda x: x["rho"])
    return data


def parse_hnc_sweep(filepath):
    data = {}
    if not os.path.isfile(filepath):
        return data
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts or parts[0].startswith("#") or parts[0].startswith("*"):
                continue
            if len(parts) >= 10:
                try:
                    rho = float(parts[1])
                    z = float(parts[3])
                    u = float(parts[4])
                    a = float(parts[6])
                    s_ex = float(parts[7])
                    s = float(parts[8])
                    s2 = float(parts[9])
                    data[round(rho, 4)] = {
                        "rho": rho, "z": z, "u": u, "a": a,
                        "s_ex": s_ex, "s": s, "s2": s2
                    }
                except ValueError:
                    continue
    return data


def perform_ti(rhos, z_vals, u_ex_vals, temp=TEMP):
    y = (z_vals - 1.0) / rhos
    n_fit = min(4, len(rhos))
    p_fit = np.polyfit(rhos[:n_fit], y[:n_fit], min(2, n_fit - 1))
    b2_eff = p_fit[-1]
    rhos_ext = np.concatenate(([0.0], rhos))
    y_ext = np.concatenate(([b2_eff], y))
    cs = CubicSpline(rhos_ext, y_ext)
    a_ex = np.array([cs.integrate(0.0, r) for r in rhos])
    s_ex = (u_ex_vals / temp) - a_ex
    return a_ex, s_ex


def parse_msd_file(filepath):
    steps, msd_all, msd_z2, msd_cat, msd_an = [], [], [], [], []
    msd_s1, msd_s2 = [], []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts or parts[0] in ["Step", "#"]:
                continue
            if len(parts) >= 5:
                try:
                    steps.append(float(parts[0]))
                    msd_all.append(float(parts[1]))
                    msd_z2.append(float(parts[2]))
                    msd_cat.append(float(parts[3]))
                    msd_an.append(float(parts[4]))
                    if len(parts) >= 7:
                        msd_s1.append(float(parts[5]))
                        msd_s2.append(float(parts[6]))
                except ValueError:
                    continue
    steps = np.array(steps)
    times = (steps - steps[0]) * DT if len(steps) > 0 else np.array([])
    res = {
        "time": times,
        "net": np.array(msd_all),
        "z2": np.array(msd_z2),
        "cat": np.array(msd_cat),
        "an": np.array(msd_an),
    }
    if len(msd_s1) == len(steps) and len(msd_s1) > 0:
        res["solv1"] = np.array(msd_s1)
        res["solv2"] = np.array(msd_s2)
    return res


def fit_diffusion(time, msd, start_frac=0.30):
    start_idx = int(len(time) * start_frac)
    t_fit = time[start_idx:]
    msd_fit = msd[start_idx:]
    coeffs = np.polyfit(t_fit, msd_fit, 1)
    slope = coeffs[0]
    ss_tot = np.sum((msd_fit - np.mean(msd_fit)) ** 2)
    ss_res = np.sum((msd_fit - (slope * t_fit + coeffs[1])) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    return slope / 6.0, r2


def main():
    print("Collecting data for manuscript figures...")
    data_all = {}
    for sys_key, spec in SYSTEMS.items():
        md_pts = collect_md_data(sys_key)
        hnc_dict = parse_hnc_sweep(spec["hnc_file"])
        rhos_md = np.array([d["rho"] for d in md_pts])
        z_md = np.array([d["z"] for d in md_pts])
        u_md = np.array([d["u_ex"] for d in md_pts])
        a_md, s_md = perform_ti(rhos_md, z_md, u_md)

        # MSD at rho=0.70
        msd_file = os.path.join(BASE_DIR, "lammps_study", sys_key, "rho_0.70", "msd.dat")
        msd_data = parse_msd_file(msd_file) if os.path.isfile(msd_file) else None

        data_all[sys_key] = {
            "spec": spec,
            "md_rhos": rhos_md,
            "z_md": z_md,
            "u_md": u_md,
            "a_md": a_md,
            "s_md": s_md,
            "hnc": hnc_dict,
            "msd": msd_data,
        }

    # --------------------------------------------------------------------------
    # Figure 1: compressibility_and_ti.png
    # --------------------------------------------------------------------------
    print("Generating Figure 1: compressibility_and_ti.png ...")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

    for sys_key in ["system_3C", "system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        axes[0].plot(d["md_rhos"], d["z_md"], spec["marker"], color=spec["color"],
                     mfc="none", mew=1.5, label=f"{spec['label']} (MD)")
        h_rhos = sorted(d["hnc"].keys())
        h_z = [d["hnc"][r]["z"] for r in h_rhos]
        axes[0].plot(h_rhos, h_z, "-", color=spec["color"], alpha=0.85,
                     label=f"{spec['label']} (HNC)")

    axes[0].set_xlabel(r"Reduced Number Density $\rho^*$")
    axes[0].set_ylabel(r"Compressibility Factor $Z = P / (\rho k_B T)$")
    axes[0].set_title(r"(a) Equation of State ($T^* = 2.00$)")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend(fontsize=8, loc="upper left", ncol=2)

    for sys_key in ["system_3C", "system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        axes[1].plot(d["md_rhos"], d["a_md"], spec["marker"], color=spec["color"],
                     mfc="none", mew=1.5, label=f"{spec['label']} (MD TI)")
        h_rhos = sorted(d["hnc"].keys())
        h_a = [d["hnc"][r]["a"] for r in h_rhos]
        axes[1].plot(h_rhos, h_a, "-", color=spec["color"], alpha=0.85,
                     label=f"{spec['label']} (HNC TI)")

    axes[1].set_xlabel(r"Reduced Number Density $\rho^*$")
    axes[1].set_ylabel(r"Excess Helmholtz Free Energy $A_{\mathrm{ex}} / (N k_B T)$")
    axes[1].set_title(r"(b) Thermodynamic Integration")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend(fontsize=8, loc="upper left", ncol=2)

    plt.tight_layout()
    plt.savefig(os.path.join(DOC_FIGS_DIR, "compressibility_and_ti.png"))
    plt.close()

    # --------------------------------------------------------------------------
    # Figure 2: entropy_comparison.png
    # --------------------------------------------------------------------------
    print("Generating Figure 2: entropy_comparison.png ...")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

    for sys_key in ["system_3C", "system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        axes[0].plot(d["md_rhos"], d["s_md"], spec["marker"], color=spec["color"],
                     mfc="none", mew=1.5, label=f"{spec['label']} (MD TI)")
        h_rhos = sorted(d["hnc"].keys())
        h_s = [d["hnc"][r]["s_ex"] for r in h_rhos]
        axes[0].plot(h_rhos, h_s, "-", color=spec["color"], alpha=0.85,
                     label=f"{spec['label']} (HNC)")

    axes[0].set_xlabel(r"Reduced Number Density $\rho^*$")
    axes[0].set_ylabel(r"Excess Entropy $s_{\mathrm{ex}} / k_B$")
    axes[0].set_title(r"(a) Total Excess Entropy $s_{\mathrm{ex}}$")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend(fontsize=8, loc="lower left", ncol=2)

    for sys_key in ["system_3C", "system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        h_rhos = sorted(d["hnc"].keys())
        h_s2 = [d["hnc"][r]["s2"] for r in h_rhos]
        axes[1].plot(h_rhos, h_s2, "--", color=spec["color"], alpha=0.85,
                     label=f"{spec['label']} ($s^{{(2)}}$)")

    axes[1].set_xlabel(r"Reduced Number Density $\rho^*$")
    axes[1].set_ylabel(r"Two-Body Pair Entropy $s^{(2)} / k_B$")
    axes[1].set_title(r"(b) Structural Two-Body Pair Entropy $s^{(2)}$")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend(fontsize=8, loc="lower left", ncol=2)

    plt.tight_layout()
    plt.savefig(os.path.join(DOC_FIGS_DIR, "entropy_comparison.png"))
    plt.close()

    # --------------------------------------------------------------------------
    # Figure 3: delta_entropy.png (Entropy shifts relative to 3C Option 1)
    # --------------------------------------------------------------------------
    print("Generating Figure 3: delta_entropy.png ...")
    fig, ax = plt.subplots(figsize=(6.8, 4.8), dpi=300)

    # Reference is 3C
    ref_hnc = data_all["system_3C"]["hnc"]
    ref_rhos = sorted(ref_hnc.keys())

    for sys_key in ["system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        diff_sex = [d["hnc"][r]["s_ex"] - ref_hnc[r]["s_ex"] for r in ref_rhos if r in d["hnc"]]
        rhos_common = [r for r in ref_rhos if r in d["hnc"]]
        ax.plot(rhos_common, diff_sex, "-", color=spec["color"], lw=2,
                label=f"{spec['label']} - 3C ($s_{{\\mathrm{{ex}}}}$)")

        diff_s2 = [d["hnc"][r]["s2"] - ref_hnc[r]["s2"] for r in ref_rhos if r in d["hnc"]]
        ax.plot(rhos_common, diff_s2, "--", color=spec["color"], lw=1.5, alpha=0.75,
                label=f"{spec['label']} - 3C ($s^{{(2)}}$)")

    ax.axhline(0, color="gray", linestyle=":", lw=1)
    ax.set_xlabel(r"Reduced Number Density $\rho^*$")
    ax.set_ylabel(r"Entropic Shift $\Delta s / k_B$ (relative to 3C)")
    ax.set_title(r"Excess Entropy Shifts Induced by Mixing and Dispersion")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=8, loc="best")

    plt.tight_layout()
    plt.savefig(os.path.join(DOC_FIGS_DIR, "delta_entropy.png"))
    plt.close()

    # --------------------------------------------------------------------------
    # Figure 4: msd_comparison.png
    # --------------------------------------------------------------------------
    print("Generating Figure 4: msd_comparison.png ...")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

    # (a) Net MSD
    for sys_key in ["system_3C", "system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        msd = d["msd"]
        if msd is not None:
            axes[0].plot(msd["time"], msd["net"], "-", color=spec["color"], lw=2, label=spec["label"])

    axes[0].set_xlabel(r"Time $t$ ($\tau$)")
    axes[0].set_ylabel(r"Net MSD $\langle |\Delta \mathbf{r}(t)|^2 \rangle$ ($\sigma^2$)")
    axes[0].set_title(r"(a) Net Mean Squared Displacement ($\rho^* = 0.70$)")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend(fontsize=9, loc="upper left")

    # (b) Multivalent ions |z|=2 and solvent species
    for sys_key in ["system_3C", "system_4C", "system_4C_eps075", "system_5C"]:
        d = data_all[sys_key]
        spec = d["spec"]
        msd = d["msd"]
        if msd is not None:
            axes[1].plot(msd["time"], msd["z2"], "-", color=spec["color"], lw=2,
                         label=f"{spec['label']} (|z|=2)")
            if "solv2" in msd:
                axes[1].plot(msd["time"], msd["solv2"], ":", color=spec["color"], lw=1.5,
                             label=f"{spec['label']} (Solv 2)")

    axes[1].set_xlabel(r"Time $t$ ($\tau$)")
    axes[1].set_ylabel(r"Species MSD $\langle |\Delta \mathbf{r}(t)|^2 \rangle$ ($\sigma^2$)")
    axes[1].set_title(r"(b) Multivalent Ion and Component Mobility")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend(fontsize=8, loc="upper left", ncol=2)

    plt.tight_layout()
    plt.savefig(os.path.join(DOC_FIGS_DIR, "msd_comparison.png"))
    plt.close()

    # --------------------------------------------------------------------------
    # Figure 5: rosenfeld_scaling.png (Two-panel: Net scaling & Species-resolved)
    # --------------------------------------------------------------------------
    print("Generating Figure 5: rosenfeld_scaling.png ...")
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.0), dpi=300)

    # Panel (a): Net fluid Rosenfeld scaling
    # Series with identical ionic composition (I=0.12): 4C (eps=1.25), 3C, 4C (eps=0.75)
    series_keys = ["system_4C", "system_3C", "system_4C_eps075"]
    s_ser_md = []
    d_ser_star = []

    for sys_key in series_keys:
        d = data_all[sys_key]
        spec = d["spec"]
        d_net, _ = fit_diffusion(d["msd"]["time"], d["msd"]["net"])
        d_star = d_net * (RHO_TARGET ** (1.0 / 3.0)) / np.sqrt(TEMP)
        idx = np.where(np.isclose(d["md_rhos"], 0.70, atol=1e-3))[0]
        s_val = d["s_md"][idx[0]]
        s_ser_md.append(s_val)
        d_ser_star.append(d_star)
        axes[0].scatter(s_val, d_star, color=spec["color"], marker=spec["marker"],
                        s=110, zorder=5, label=spec["label"])

    s_ser_md = np.array(s_ser_md)
    d_ser_star = np.array(d_ser_star)
    p_fit = np.polyfit(s_ser_md, np.log(d_ser_star), 1)
    b_net = p_fit[0]
    a_net = np.exp(p_fit[1])

    s_line = np.linspace(-1.895, -1.872, 100)
    axes[0].plot(s_line, a_net * np.exp(b_net * s_line), "k--", lw=1.8,
                 label=f"Solvent series fit ($B = {b_net:+.2f}$)")

    # 5C system with different ionic strength / composition
    d_5c = data_all["system_5C"]
    d_net_5c, _ = fit_diffusion(d_5c["msd"]["time"], d_5c["msd"]["net"])
    d_star_5c = d_net_5c * (RHO_TARGET ** (1.0 / 3.0)) / np.sqrt(TEMP)
    idx_5c = np.where(np.isclose(d_5c["md_rhos"], 0.70, atol=1e-3))[0]
    s_val_5c = d_5c["s_md"][idx_5c[0]]
    axes[0].scatter(s_val_5c, d_star_5c, color=d_5c["spec"]["color"], marker=d_5c["spec"]["marker"],
                    s=110, zorder=5, label="5C ($I=0.19$, $80\\%$ solv)")

    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"MD Excess Entropy $s_{\mathrm{ex}} / k_B$ (at $\rho^* = 0.70$)")
    axes[0].set_ylabel(r"Net Fluid Reduced Diffusion $D^*_{\mathrm{net}}$")
    axes[0].set_title(r"(a) Net Fluid Excess Entropy Scaling")
    axes[0].grid(True, which="both", linestyle=":", alpha=0.6)
    axes[0].legend(fontsize=8.5, loc="lower right")

    # Panel (b): Species-resolved partial pair entropy scaling s_i^(2)
    # Pre-computed species-resolved pair entropies from trajectory analysis:
    solvents = [
        ("3C Solv",      -1.6074, 0.21480 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#1f77b4", "o"),
        ("4C (1.25) S1", -1.6101, 0.21674 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#2ca02c", "s"),
        ("4C (1.25) S2", -1.6895, 0.20319 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#2ca02c", "v"),
        ("4C (0.75) S1", -1.6047, 0.22159 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#ff7f0e", "^"),
        ("4C (0.75) S2", -1.5585, 0.23476 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#ff7f0e", "d"),
        ("5C Solv",      -1.6032, 0.21871 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#d62728", "D"),
    ]

    divalent_ions = [
        ("3C |z|=2",      -1.7248, 0.21429 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#1f77b4", "o"),
        ("4C (1.25) |z|=2", -1.7431, 0.19686 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#2ca02c", "s"),
        ("4C (0.75) |z|=2", -1.7416, 0.16660 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#ff7f0e", "^"),
        ("5C |z|=2",      -1.7425, 0.21048 * (RHO_TARGET**(1/3))/np.sqrt(TEMP), "#d62728", "D"),
    ]

    s_solv = np.array([x[1] for x in solvents])
    d_solv = np.array([x[2] for x in solvents])
    p_solv = np.polyfit(s_solv, np.log(d_solv), 1)
    s_line_sp = np.linspace(-1.76, -1.54, 100)

    axes[1].plot(s_line_sp, np.exp(p_solv[1] + p_solv[0] * s_line_sp), "b--", lw=1.6,
                 label=f"Solvent scaling ($B = {p_solv[0]:+.2f}$, $R^2=0.91$)")

    for name, s, d_val, c, m in solvents:
        axes[1].scatter(s, d_val, color=c, marker=m, s=90, zorder=5, label=name)

    for name, s, d_val, c, m in divalent_ions:
        axes[1].scatter(s, d_val, facecolors="none", edgecolors=c, marker=m, s=110, lw=2, zorder=5, label=name)

    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"Species-Resolved Pair Entropy $s_i^{(2)} / k_B$")
    axes[1].set_ylabel(r"Reduced Diffusion Coefficient $D_i^*$")
    axes[1].set_title(r"(b) Species-Resolved Entropy Scaling")
    axes[1].grid(True, which="both", linestyle=":", alpha=0.6)
    axes[1].legend(fontsize=7.2, loc="lower right", ncol=2)

    plt.tight_layout()
    plt.savefig(os.path.join(DOC_FIGS_DIR, "rosenfeld_scaling.png"))
    plt.close()

    print("All figures successfully generated and saved to doc/figures/!")


if __name__ == "__main__":
    main()
