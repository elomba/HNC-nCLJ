#!/usr/bin/env python3
"""
Script: thermodynamic_integration.py
Description: Performs thermodynamic integration along the density isotherm:
             A_ex(rho)/(N*k_B*T) = \int_0^rho [(Z(rho') - 1)/rho'] d rho'
             and compares the excess entropy s_ex and two-body entropy s^(2)
             between the 3-Component and 5-Component electrolyte systems.
"""

import os
import sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def parse_thermohncr(filepath):
    """Parses Thermohncr.dat file into a structured dictionary of arrays."""
    data_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 12:
                data_lines.append([float(p) for p in parts])

    if not data_lines:
        raise ValueError(f"No numerical data found in {filepath}")

    arr = np.array(data_lines)
    # Sort by density ascending
    arr = arr[arr[:, 1].argsort()]

    # Column mapping based on format:
    # 0: x(1), 1: rho, 2: beta*, 3: z, 4: u/nkt, 5: u_c/nkt, 6: a/nkt,
    # 7: s_ex/nk, 8: s/nk, 9: s2/nk, 10: z(a), 11: xc
    res = {
        'x1': arr[:, 0],
        'rho': arr[:, 1],
        'beta': arr[:, 2],
        'z_vir': arr[:, 3],
        'u_ex': arr[:, 4],
        'u_c': arr[:, 5],
        'a_hnc': arr[:, 6],
        's_ex': arr[:, 7],
        's_tot': arr[:, 8],
        's2': arr[:, 9],
        'z_a': arr[:, 10],
        'xc': arr[:, 11],
    }
    return res

def perform_thermodynamic_integration(rho, z_vir):
    """
    Computes excess free energy by integrating (Z - 1) / rho from 0 to rho.
    Extrapolates (Z - 1)/rho to rho -> 0 using polynomial fit.
    """
    y = (z_vir - 1.0) / rho
    # Fit quadratic polynomial to low-density points to extrapolate to rho=0 (second virial coefficient)
    p_fit = np.polyfit(rho[:4], y[:4], 2)
    b2_eff = p_fit[-1]

    # Prepend rho=0 with y(0) = B2
    rho_full = np.concatenate(([0.0], rho))
    y_full = np.concatenate(([b2_eff], y))

    # Cubic spline integration
    cs = CubicSpline(rho_full, y_full)
    a_ti = np.array([cs.integrate(0.0, r) for r in rho])
    return a_ti, b2_eff

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.abspath(os.path.join(script_dir, '..', 'results'))
    os.makedirs(results_dir, exist_ok=True)

    file_3c = os.path.join(results_dir, 'Thermo_3C_density_sweep.dat')
    file_5c = os.path.join(results_dir, 'Thermo_5C_density_sweep.dat')

    if not os.path.isfile(file_3c) or not os.path.isfile(file_5c):
        print(f"Error: Missing input files in {results_dir}")
        print(f"Looked for: {file_3c} and {file_5c}")
        sys.exit(1)

    d3 = parse_thermohncr(file_3c)
    d5 = parse_thermohncr(file_5c)

    # Perform thermodynamic integration
    a_ti_3c, b2_3c = perform_thermodynamic_integration(d3['rho'], d3['z_vir'])
    a_ti_5c, b2_5c = perform_thermodynamic_integration(d5['rho'], d5['z_vir'])

    # Entropy derived from thermodynamic integration
    s_ex_ti_3c = d3['u_ex'] - a_ti_3c
    s_ex_ti_5c = d5['u_ex'] - a_ti_5c

    # Differences (5C - 3C)
    delta_s_ex = d5['s_ex'] - d3['s_ex']
    delta_s2 = d5['s2'] - d3['s2']
    delta_u_ex = d5['u_ex'] - d3['u_ex']
    delta_u_c = d5['u_c'] - d3['u_c']
    delta_z = d5['z_vir'] - d3['z_vir']

    # Text summary output
    summary_path = os.path.join(results_dir, 'comparison_summary.txt')
    with open(summary_path, 'w') as f:
        header = (
            "====================================================================================================\n"
            "   COMPARATIVE ENTROPY & THERMODYNAMIC ANALYSIS: 3-COMPONENT vs 5-COMPONENT CHARGED FLUIDS\n"
            "   Conditions: beta* = 0.50 (T* = 2.00), gamma = 1.0, Solvent x_1 = 0.80, Total Ions = 0.20\n"
            "====================================================================================================\n\n"
            f"Extrapolated Effective Second Virial Coefficient B2:\n"
            f"  System 3C (Solvent + 2:2 Salt): B2* = {b2_3c:+.4f}\n"
            f"  System 5C (Solvent + 1:1 & 2:2 Salt): B2* = {b2_5c:+.4f}\n\n"
            "----------------------------------------------------------------------------------------------------\n"
            " rho   |   Z_3C    Z_5C  |  s_ex_3C s_ex_5C Delta_sex |   s2_3C   s2_5C  Delta_s2 |   s2/sex(3C, 5C)  \n"
            "----------------------------------------------------------------------------------------------------\n"
        )
        f.write(header)
        print(header, end='')

        for i, r in enumerate(d3['rho']):
            r3_ratio = d3['s2'][i] / d3['s_ex'][i] if abs(d3['s_ex'][i]) > 1e-4 else 1.0
            r5_ratio = d5['s2'][i] / d5['s_ex'][i] if abs(d5['s_ex'][i]) > 1e-4 else 1.0
            line = (
                f"{r:5.2f} | {d3['z_vir'][i]:7.4f} {d5['z_vir'][i]:7.4f} | "
                f"{d3['s_ex'][i]:8.4f} {d5['s_ex'][i]:8.4f} {delta_s_ex[i]:+9.4f} | "
                f"{d3['s2'][i]:7.4f} {d5['s2'][i]:7.4f} {delta_s2[i]:+9.4f} | "
                f"{r3_ratio*100:6.1f}%  {r5_ratio*100:6.1f}%\n"
            )
            f.write(line)
            print(line, end='')

        f.write("----------------------------------------------------------------------------------------------------\n")
        print("----------------------------------------------------------------------------------------------------")

        ti_check_header = (
            "\nConsistency Check: Direct HNC k-space Free Energy (a_hnc) vs Thermodynamic Integration (a_ti):\n"
            "----------------------------------------------------------------------------------------------------\n"
            " rho   |  a_hnc(3C)   a_ti(3C)   Diff(3C)  |  a_hnc(5C)   a_ti(5C)   Diff(5C)  \n"
            "----------------------------------------------------------------------------------------------------\n"
        )
        f.write(ti_check_header)
        print(ti_check_header, end='')

        for i, r in enumerate(d3['rho']):
            diff3 = d3['a_hnc'][i] - a_ti_3c[i]
            diff5 = d5['a_hnc'][i] - a_ti_5c[i]
            line = (
                f"{r:5.2f} | {d3['a_hnc'][i]:10.4f} {a_ti_3c[i]:10.4f} {diff3:+10.4f} | "
                f"{d5['a_hnc'][i]:10.4f} {a_ti_5c[i]:10.4f} {diff5:+10.4f}\n"
            )
            f.write(line)
            print(line, end='')
        f.write("----------------------------------------------------------------------------------------------------\n")
        print("----------------------------------------------------------------------------------------------------")

    # CSV output
    csv_path = os.path.join(results_dir, 'comparison_data.csv')
    csv_data = np.column_stack([
        d3['rho'],
        d3['z_vir'], d5['z_vir'],
        d3['u_ex'], d5['u_ex'],
        d3['u_c'], d5['u_c'],
        d3['a_hnc'], d5['a_hnc'],
        a_ti_3c, a_ti_5c,
        d3['s_ex'], d5['s_ex'], delta_s_ex,
        d3['s2'], d5['s2'], delta_s2,
        d3['xc'], d5['xc']
    ])
    csv_header = (
        "rho,z_3c,z_5c,u_ex_3c,u_ex_5c,u_c_3c,u_c_5c,a_hnc_3c,a_hnc_5c,"
        "a_ti_3c,a_ti_5c,s_ex_3c,s_ex_5c,delta_s_ex,s2_3c,s2_5c,delta_s2,xc_3c,xc_5c"
    )
    np.savetxt(csv_path, csv_data, delimiter=',', header=csv_header, comments='', fmt='%.6f')
    print(f"\nSaved CSV data: {csv_path}")

    # Generate publication-quality plots
    if HAS_MATPLOTLIB:
        plt.rcParams.update({
            'font.size': 11,
            'axes.labelsize': 12,
            'legend.fontsize': 10,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'figure.dpi': 300
        })

        # Figure 1: Excess entropy s_ex and pair entropy s^(2)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(d3['rho'], d3['s_ex'], 'b-o', label=r'$s_{\mathrm{ex}}$ (3C: pure 2:2)', lw=2, ms=5)
        ax.plot(d5['rho'], d5['s_ex'], 'r-s', label=r'$s_{\mathrm{ex}}$ (5C: 1:1 + 2:2)', lw=2, ms=5)
        ax.plot(d3['rho'], d3['s2'], 'b--', label=r'$s^{(2)}$ (3C)', lw=1.5)
        ax.plot(d5['rho'], d5['s2'], 'r--', label=r'$s^{(2)}$ (5C)', lw=1.5)
        ax.set_xlabel(r'Total Reduced Density $\rho^*$')
        ax.set_ylabel(r'Entropy per Particle $[k_B]$')
        ax.set_title(r'Excess and Two-Body Entropy vs Density ($\beta^* = 0.50, \gamma = 1.0$)')
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True)
        fig.tight_layout()
        p1 = os.path.join(results_dir, 'entropy_comparison.png')
        fig.savefig(p1)
        plt.close(fig)

        # Figure 2: Delta s_ex and Delta s2 (5C - 3C)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(d3['rho'], delta_s_ex, 'g-o', label=r'$\Delta s_{\mathrm{ex}} = s_{\mathrm{ex}}^{\mathrm{5C}} - s_{\mathrm{ex}}^{\mathrm{3C}}$', lw=2, ms=5)
        ax.plot(d3['rho'], delta_s2, 'm--s', label=r'$\Delta s^{(2)} = s^{(2),\mathrm{5C}} - s^{(2),\mathrm{3C}}$', lw=2, ms=5)
        ax.axhline(0, color='gray', linestyle=':', lw=1)
        ax.set_xlabel(r'Total Reduced Density $\rho^*$')
        ax.set_ylabel(r'Entropy Difference $[k_B]$')
        ax.set_title(r'Entropy Difference (5C Mixture vs 3C Pure 2:2 Electrolyte)')
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True)
        fig.tight_layout()
        p2 = os.path.join(results_dir, 'delta_entropy.png')
        fig.savefig(p2)
        plt.close(fig)

        # Figure 3: Compressibility and Thermodynamic Integration
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.plot(d3['rho'], d3['z_vir'], 'b-o', label='3C (virial)', lw=1.5, ms=4)
        ax1.plot(d5['rho'], d5['z_vir'], 'r-s', label='5C (virial)', lw=1.5, ms=4)
        ax1.plot(d3['rho'], d3['z_a'], 'b--', label=r'3C (from $a_{\mathrm{ex}}$)', lw=1.2)
        ax1.plot(d5['rho'], d5['z_a'], 'r--', label=r'5C (from $a_{\mathrm{ex}}$)', lw=1.2)
        ax1.set_xlabel(r'$\rho^*$')
        ax1.set_ylabel(r'Compressibility Factor $Z = P / (\rho k_B T)$')
        ax1.set_title(r'Compressibility Factor $Z(\rho^*)$')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(d3['rho'], d3['a_hnc'], 'b-o', label=r'3C (direct $a_{\mathrm{HNC}}$)', lw=1.5, ms=4)
        ax2.plot(d3['rho'], a_ti_3c, 'b--', label=r'3C (integrated $A_{\mathrm{TI}}$)', lw=1.5)
        ax2.plot(d5['rho'], d5['a_hnc'], 'r-s', label=r'5C (direct $a_{\mathrm{HNC}}$)', lw=1.5, ms=4)
        ax2.plot(d5['rho'], a_ti_5c, 'r--', label=r'5C (integrated $A_{\mathrm{TI}}$)', lw=1.5)
        ax2.set_xlabel(r'$\rho^*$')
        ax2.set_ylabel(r'Excess Free Energy $A_{\mathrm{ex}} / (N k_B T)$')
        ax2.set_title(r'Free Energy: Direct vs Integrated')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        fig.tight_layout()
        p3 = os.path.join(results_dir, 'compressibility_and_ti.png')
        fig.savefig(p3)
        plt.close(fig)

        print(f"Generated plots:\n  {p1}\n  {p2}\n  {p3}")

if __name__ == '__main__':
    main()
