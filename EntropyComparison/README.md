# Comparative Study: Excess Entropy and Structural Ordering in 3-Component vs. 5-Component Charged LJ Fluids

This project investigates the configurational and excess entropy differences between two multicomponent ionic fluids with identical solvent concentration and packing, subjected to thermodynamic integration along a supercritical density path down to $\rho^* = 0.01$.

---

## 1. Physical Model & Systems Compared

Both systems share identical packing fractions ($\eta \approx 0.3665$ at $\rho^* = 0.70$), identical solvent mole fraction ($x_{\text{solvent}} = 0.80$), identical total ionic solute fraction ($x_{\text{ions}} = 0.20$), and identical short-range Lennard-Jones parameters ($\sigma_{ij} = 1.0$, $\epsilon_{ij} = 1.0$, $r_{\text{cut},ij} = 2.5$ for all pairs).

| Parameter | System 3C (Solvent + 2:2 Electrolyte) | System 5C (Solvent + 1:1 \& 2:2 Mixed Electrolyte) |
|:---|:---|:---|
| **Species Count** | $NSP = 3$ | $NSP = 5$ |
| **Species 1 (Solvent)** | $q_1 = 0.0$, $x_1 = 0.80$ | $q_1 = 0.0$, $x_1 = 0.80$ |
| **Species 2** | Monovalent-free, Divalent Cation $q_2 = +2.0$, $x_2 = 0.10$ | Monovalent Cation $q_2 = +1.0$, $x_2 = 0.05$ |
| **Species 3** | Divalent Anion $q_3 = -2.0$, $x_3 = 0.10$ | Divalent Cation $q_3 = +2.0$, $x_3 = 0.04$ |
| **Species 4** | --- | Monovalent Anion $q_4 = -1.0$, $x_4 = 0.09$ |
| **Species 5** | --- | Divalent Anion $q_5 = -2.0$, $x_5 = 0.02$ |
| **Total Solvent $x_{\text{solv}}$** | $0.80$ | $0.80$ |
| **Total Ions $x_{\text{ions}}$** | $0.20$ | $0.20$ |
| **Ionic Strength** $I = \frac{1}{2}\sum x_i q_i^2$ | **$0.40$** | **$0.19$** |
| **Packing $\eta$ at $\rho^* = 0.70$** | $0.3665$ | $0.3665$ |
| **Coupling $\gamma$** | $1.0$ | $1.0$ |

---

## 2. Selection of Supercritical Temperature ($T^* = 2.0$, $\beta^* = 0.50$)

For pure Lennard-Jones fluids, the vapor-liquid critical temperature is $T_c^* \approx 1.31$ ($\beta_c^* \approx 0.76$). 

By fixing the temperature at:
$$T^* = 2.00 \quad (\beta^* = 0.50)$$
the fluid is maintained safely above the vapor-liquid critical point across the entire density range $\rho^* \in [0.01, 0.70]$. Consequently:
1. The isothermal compressibility $\chi_T^{-1}(\rho^*) = \frac{1}{k_B T}\left(\frac{\partial P}{\partial \rho}\right) > 0$ remains strictly positive and bounded away from zero along the entire path.
2. The system never encounters a spinodal instability or two-phase boundary.
3. The HNC integral equations converge stably and rapidly at every density point via Kinoshita--Harada continuation.

---

## 3. Thermodynamic Integration

Along the chosen isotherm ($T^* = 2.0$, $\beta^* = 0.50$), the excess Helmholtz free energy is obtained by integrating the virial equation of state from the ideal gas limit ($\rho \to 0$):
$$\frac{A_{\text{ex}}(\rho)}{N k_B T} = \int_0^{\rho} \frac{Z(\rho') - 1}{\rho'} d\rho'$$

where $Z(\rho) = P / (\rho k_B T)$ is the compressibility factor.

The excess entropy is then evaluated both directly from HNC:
$$\frac{s_{\text{ex}}^{\text{HNC}}}{N k_B} = \frac{u_{\text{ex}} - a_{\text{ex}}}{N k_B T}$$
and via thermodynamic integration:
$$\frac{s_{\text{ex}}^{\text{TI}}}{N k_B} = \frac{u_{\text{ex}}}{N k_B T} - \frac{A_{\text{ex}}^{\text{TI}}}{N k_B T}$$

---

## 4. Project Directory Structure

```
.
├── system_3C/
│   ├── input_init_rho07.hncr       # Thermal sweep to reach rho=0.70, beta*=0.50
│   └── input_density_sweep.hncr    # Continuation density sweep (0.70 -> 0.01)
├── system_5C/
│   ├── input_init_rho07.hncr       # Thermal sweep to reach rho=0.70, beta*=0.50
│   └── input_density_sweep.hncr    # Continuation density sweep (0.70 -> 0.01)
├── scripts/
│   ├── run_study.sh                # End-to-end execution pipeline
│   └── thermodynamic_integration.py# TI integration, post-processing, and plotting
└── results/
    ├── Thermo_3C_density_sweep.dat # Density sweep outputs for 3C
    ├── Thermo_5C_density_sweep.dat # Density sweep outputs for 5C
    ├── comparison_summary.txt      # Text summary table
    ├── comparison_data.csv         # Complete tabular dataset
    ├── entropy_comparison.png      # Plot: s_ex and s^(2) vs density
    ├── delta_entropy.png           # Plot: Delta s_ex and Delta s^(2) vs density
    └── compressibility_and_ti.png  # Plot: Z(rho) and free energy consistency
```

---

## 5. How to Run the Complete Study

From the `EntropyComparison/scripts/` directory:
```bash
cd EntropyComparison/scripts
./run_study.sh
```
This automatically compiles the latest binary if needed, runs the initialization and density continuation sweeps for both systems, carries out the numerical thermodynamic integration, and writes plots and summary tables to `results/`.
