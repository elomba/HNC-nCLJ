#!/usr/bin/env python3
"""
Calculate HNC partial molar excess properties for the corrected System 5C:
x = [0.84, 0.05, 0.03, 0.05, 0.03]
"""

import os, sys, subprocess
import numpy as np

BIN = "/home/e.lomba/EntropicFLuids/bin/hncr_nclj"
WORK_DIR = "/home/e.lomba/EntropicFLuids/EntropyComparison/system_5C_partial"
os.makedirs(WORK_DIR, exist_ok=True)

# 1. First run beta sweep at rho = 0.7000: beta in [0.4950, 0.5000, 0.5050]
input_beta = f"""5-Component System (Solvent + 1:1 and 2:2 Electrolytes): Beta Sweep
4096 0.005 5.0d-7 0.5 10000.0 0.05
100 5.0d-7
5
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
0.0d0 1.0d0 2.0d0 -1.0d0 -2.0d0
1 3
0.7000
0.4950 0.5000 0.5050
.f.
0.84 0.05 0.03
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
1.0 1
"""

# Copy converged solin from EntropyComparison/system_5C
os.system(f"cp /home/e.lomba/EntropicFLuids/EntropyComparison/system_5C/lsol.hncr {WORK_DIR}/solin.hncr")

with open(f"{WORK_DIR}/input_beta.hncr", "w") as f:
    f.write(input_beta)

print("Running beta sweep for System 5C...")
subprocess.run([BIN, "input_beta.hncr"], cwd=WORK_DIR, check=True)
os.system(f"cp {WORK_DIR}/Thermohncr.dat {WORK_DIR}/Thermo_beta.dat")

# 2. Now run rho sweep at beta = 0.5000: rho in [0.6950, 0.7000, 0.7050]
input_rho = f"""5-Component System (Solvent + 1:1 and 2:2 Electrolytes): Rho Sweep
4096 0.005 5.0d-7 0.5 10000.0 0.05
100 5.0d-7
5
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0 1.0
0.0d0 1.0d0 2.0d0 -1.0d0 -2.0d0
3 1
0.6950 0.7000 0.7050
0.5000
.f.
0.84 0.05 0.03
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
2.5 2.5 2.5 2.5 2.5
1.0 1
"""

with open(f"{WORK_DIR}/input_rho.hncr", "w") as f:
    f.write(input_rho)

print("Running rho sweep for System 5C...")
subprocess.run([BIN, "input_rho.hncr"], cwd=WORK_DIR, check=True)
os.system(f"cp {WORK_DIR}/Thermohncr.dat {WORK_DIR}/Thermo_rho.dat")

# 3. Parse and compute derivatives
def parse_thermo(filename):
    data = []
    with open(filename) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 17 and parts[0] != "#" and parts[0] != "x(1)":
                try:
                    vals = [float(p) for p in parts]
                    data.append(vals)
                except:
                    pass
    return np.array(data)

d_beta = parse_thermo(f"{WORK_DIR}/Thermo_beta.dat")
d_rho = parse_thermo(f"{WORK_DIR}/Thermo_rho.dat")

# Columns in Thermohncr.dat:
# 0: x(1)
# 1: rho
# 2: beta*
# 3: z
# 4: u/nkt
# 5: u_c/nkt
# 6: a/nkt
# 7: s_ex/nk
# 8: s/nk
# 9: s2/nk
# 10: z(a)
# 11: xc
# 12: mu/kt(1)
# 13: mu/kt(2)
# 14: mu/kt(3)
# 15: mu/kt(4)
# 16: mu/kt(5)

# Sort beta
d_beta = d_beta[np.argsort(d_beta[:, 2])]
# Sort rho
d_rho = d_rho[np.argsort(d_rho[:, 1])]

betas = d_beta[:, 2]
rhos = d_rho[:, 1]

# P = rho * Z / beta
P_beta = d_beta[:, 1] * d_beta[:, 3] / d_beta[:, 2]
# T = 1 / beta
T_beta = 1.0 / betas
# gamma_V = (dP / dT)_V
gamma_V = (P_beta[2] - P_beta[0]) / (T_beta[0] - T_beta[2])  # note T decreases as beta increases

P_rho = d_rho[:, 1] * d_rho[:, 3] / d_rho[:, 2]
dP_drho = (P_rho[2] - P_rho[0]) / (rhos[2] - rhos[0])

print(f"\nEOS Properties at rho=0.70, T=2.00:")
print(f"  gamma_V = (dP/dT)_V = {gamma_V:.4f}")
print(f"  (dP/drho)_T = {dP_drho:.4f}")

species_names = ["Solvent", "Cation (z=+1)", "Cation (z=+2)", "Anion (z=-1)", "Anion (z=-2)"]
x_fracs = [0.84, 0.05, 0.03, 0.05, 0.03]

print(f"\n{'Species':<16} {'x_i':<6} {'mu_ex/(kT)':<12} {'V_bar':<10} {'s_i,V_ex':<12} {'S_bar_P_ex':<12} {'s2_i':<10}")
print("-" * 80)

# s2 components from RDF if needed, or from Thermohncr
s2_tot = d_beta[1, 9]

results = []
for i in range(5):
    col = 12 + i
    mu_kt_0 = d_beta[1, col]
    # y = beta * mu_ex = mu_kt
    # d(y)/d(beta)
    dy_dbeta = (d_beta[2, col] - d_beta[0, col]) / (betas[2] - betas[0])
    s_i_V = -(dy_dbeta - mu_kt_0)  # in units of k_B

    # d(mu_ex)/d(rho):
    # mu_ex = mu_kt * T = mu_kt / beta
    mu_rho = d_rho[:, col] / d_rho[:, 2]
    dmu_drho = (mu_rho[2] - mu_rho[0]) / (rhos[2] - rhos[0])

    # V_bar_i = ( (kT/rho) + dmu_ex/drho ) / (dP/drho)
    kT = 1.0 / 0.50
    V_bar = (kT / 0.70 + dmu_drho) / dP_drho

    # S_bar_P_ex = s_i_V + (gamma_V / (dP/drho)) * dmu_drho
    S_bar_P = s_i_V + (gamma_V / dP_drho) * dmu_drho

    results.append({
        "name": species_names[i],
        "x": x_fracs[i],
        "mu_kt": mu_kt_0,
        "V_bar": V_bar,
        "s_i_V": s_i_V,
        "S_bar_P": S_bar_P,
    })
    print(f"{species_names[i]:<16} {x_fracs[i]:<6.2f} {mu_kt_0:<12.4f} {V_bar:<10.3f} {s_i_V:<12.3f} {S_bar_P:<12.3f}")

V_mix = sum(r['x'] * r['V_bar'] for r in results)
print(f"\nSum x_i V_bar_i = {V_mix:.4f} (expected 1/rho = {1.0/0.70:.4f})")
