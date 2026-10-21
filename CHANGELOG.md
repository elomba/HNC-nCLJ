# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [V1.6] - 2026-10-21

### Changed
- **Symmetric 5-Component Fluid Composition & Systematic Recalibration**:
  - Corrected the mole fractions of the 5-component mixture (System 5C) to enforce strict charge symmetry and exact divalent ion concentration matching with the other mixtures ($x_{+2} = x_{-2} = 0.03$, $3\%$ Cation 2+ and $3\%$ Anion 2-):
    - Solvent: $x_1 = 0.84$ ($N_1 = 1680, q_1 = 0$)
    - Cation 1+: $x_2 = 0.05$ ($N_2 = 100, q_2 = +1$)
    - Cation 2+: $x_3 = 0.03$ ($N_3 = 60, q_3 = +2$)
    - Anion 1-: $x_4 = 0.05$ ($N_4 = 100, q_4 = -1$)
    - Anion 2-: $x_5 = 0.03$ ($N_5 = 60, q_5 = -2$)
    - Ionic strength: $I = \frac{1}{2}\sum x_i q_i^2 = 0.17$; Debye screening length: $\kappa^{-1} = 0.818\,\sigma$.
  - All four fluids (3C, 4C $\varepsilon_2=1.25$, 4C $\varepsilon_2=0.75$, 5C) now share identical divalent salt mole fractions ($3\%$ Cation 2+, $3\%$ Anion 2-).

### Added
- **Full Hypernetted Chain (HNC) Recalculation for Corrected 5C**:
  - Converged HNC solution and carried out full density continuation sweep across 15 density points ($\rho^* = 0.01 \to 0.70$) with Allnatt screening.
  - Performed HNC thermodynamic integration at $\rho^* = 0.70$: $Z_{\text{virial}} = Z_{\text{diff}} = 3.2099$, $u_{\text{ex}}/(k_B T) = -1.7403$, $A_{\text{ex}}/(N k_B T) = +0.3274$, $s_{\text{ex}} = -2.0677\,k_B$, $s^{(2)} = -1.5731\,k_B$ ($s^{(2)}/s_{\text{ex}} = 76.1\%$).
  - Evaluated species-resolved partial molar excess properties ($\mu_i^{\text{ex}}, \bar{V}_i, s_{i,V}^{\text{ex}}, \bar{S}_{i,P}^{\text{ex}}$):
    - Solvent: $\mu_1^{\text{ex}}/(k_B T) = 2.6021, \bar{V}_1 = 1.434\,\sigma^3, s_{1,V}^{\text{ex}} = -5.267\,k_B, \bar{S}_{1,P}^{\text{ex}} = -1.384\,k_B$.
    - Monovalent ions ($|z|=1$): $\mu_{|z|=1}^{\text{ex}}/(k_B T) = 2.4380, \bar{V} = 1.427\,\sigma^3, s_{V}^{\text{ex}} = -5.308\,k_B, \bar{S}_{P}^{\text{ex}} = -1.445\,k_B$.
    - Divalent ions ($|z|=2$): $\mu_{|z|=2}^{\text{ex}}/(k_B T) = 1.9421, \bar{V} = 1.408\,\sigma^3, s_{V}^{\text{ex}} = -5.442\,k_B, \bar{S}_{P}^{\text{ex}} = -1.636\,k_B$.
    - Verified strict electrostrictive volume compression: $\bar{V}_{|z|=2} (1.408\,\sigma^3) < \bar{V}_{|z|=1} (1.427\,\sigma^3) < \bar{V}_{\text{solv}} (1.434\,\sigma^3)$.
- **GPU-Accelerated MD Batch & Thermodynamic Integration for Corrected 5C**:
  - Re-generated and ran all 15 LAMMPS MD simulations across the supercritical isotherm on NVIDIA RTX PRO 4500 Blackwell GPUs.
  - Computed equation of state and thermodynamic integration: $Z_{\text{MD}} = 2.5210 \pm 0.003$, $u_{\text{ex}}/(k_B T) = -1.8355 \pm 0.001$, $A_{\text{ex}}/(N k_B T) = +0.0499$, $s_{\text{ex}}^{\text{MD}} = -1.8854\,k_B$.
  - Evaluated long-time Einstein self-diffusivities at $\rho^* = 0.70$: $D_{\text{net}} = 0.21789 \pm 0.00003$, $D_{z2} = 0.22887 \pm 0.00016$ ($D_{\text{cat2}} = 0.22854, D_{\text{an2}} = 0.22684$).
- **Two-Phase Thermodynamics (2PT) MD Simulation for Corrected 5C**:
  - Performed high-frequency velocity dumping NVT simulation (30,000 steps, 5,000 snapshots over 50.0 $\tau$).
  - Extracted 2PT dynamic decomposition: $D_{\text{GK}} = 0.2219$, $S(0) = 1.2668$, $f = 0.3592$, $s_{\text{ex}}^{\text{diff}} = -0.685\,k_B$, $s_{\text{ex}}^{\text{solid}} = -4.687\,k_B$, yielding $s_{\text{mix},V}^{\text{2PT}} = -5.372\,k_B$ (matching HNC $-5.282\,k_B$ to within $1.7\%$).
- **Manuscript & Figures Update**:
  - Regenerated all figures (`compressibility_and_ti.png`, `entropy_comparison.png`, `delta_entropy.png`, `msd_comparison.png`, `rosenfeld_scaling.png`, and `vacf_and_dos_2pt.png`).
  - Updated Tables I, II, III, and IV, and text in `doc/paper_pre.tex`.
  - Expanded Table IV to include simulated Einstein self-diffusion coefficients $D_{\text{MSD}}$ alongside Green-Kubo values $D_{\text{GK}}$ for all individual species and net fluid across all four mixtures.
  - Added in-depth discussion in Section IV.D (Item 5) clarifying the physical origin of the discrepancy between $D_{\text{MSD}}$ and $D_{\text{GK}}$ for multivalent ions: localized short-time in-cage rattling ($50\,\tau$ trajectory, $10\,\tau$ VACF window, $D_{\text{GK}} \approx 0.195$) versus long-time hydrodynamic transport requiring whole contact-pair diffusion and dissociation ($300\,\tau$ MSD trajectory, $D_{z2} = 0.214 \to 0.167$).
  - Recompiled publication PDF `doc/paper_pre.pdf` (12 pages, clean compilation).

---

## [V1.5] - 2026-10-21

### Added
- **Multi-System Two-Phase Thermodynamics (2PT) Analysis across All Fluids**:
  - Configured and executed high-frequency NVT velocity dumping MD simulations on NVIDIA RTX PRO 4500 GPUs for all remaining multicomponent mixtures: System 4C ($\varepsilon_2 = 1.25$), System 4C ($\varepsilon_2 = 0.75$), and System 5C ($I = 0.19$).
  - Developed generalized analysis framework (`lammps_study/two_pt_analysis/run_multisystem_2pt.py` and `generate_fig6.py`).
  - Evaluated species-resolved Green-Kubo diffusion coefficients $D_{\text{GK}}$, zero-frequency densities of states $S(0)$, fluidicity parameters $f$, diffusive hard-sphere excess entropies $s_{\text{ex}}^{\text{diff}}$, and solid-like vibrational excess entropies $s_{\text{ex}}^{\text{solid}}$.
  - Expanded Table IV into a comprehensive multi-system table displaying the complete species breakdown (solvents, cations, and anions) side-by-side with HNC canonical constant-volume partial excess entropies $s_{i,V}^{\text{ex}}$.
  - Demonstrated quantitative agreement between 2PT and HNC across all species in all four fluids ($0.05\%\text{--}2.4\%$ error).
  - Validated solvent dispersion modulation in 2PT: weakening solvent attraction ($\varepsilon_2 = 1.25 \to 0.75$) yields a $+0.279\,k_B$ entropic surge ($s_{2,V}^{\text{2PT}} = -5.485 \to -5.206\,k_B$), quantitatively mirroring HNC ($+0.332\,k_B$).
  - Confirmed electrostriction valence ordering in 5C: divalent ions experience significantly deeper vibrational caging and entropic confinement ($s_{i,V}^{\text{2PT}} \approx -5.44$ to $-5.46\,k_B$) than monovalent ions ($-5.39$ to $-5.42\,k_B$) and neutral solvent ($-5.38\,k_B$).

### Changed
- **Figure 6 Redesign (`doc/figures/vacf_and_dos_2pt.png`)**:
  - Upgraded Figure 6 into a publication-quality 4-panel layout:
    - **Panel (a)**: Normalized solvent velocity autocorrelation functions $C_v(t) / C_v(0)$ comparing 3C, 4C ($\varepsilon_2=1.25$), and 4C ($\varepsilon_2=0.75$).
    - **Panel (b)**: Ionic VACFs in System 5C, contrasting backscattering well depths for monovalent vs divalent ions.
    - **Panel (c)**: Solvent vibrational density of states $S(\nu)$, showing the elevated diffusive peak $S(0)$ and softened vibrational band for $\varepsilon_2 = 0.75$.
    - **Panel (d)**: Ionic vibrational density of states in System 5C, illustrating the electrostrictive blue-shift of solid-like rattling modes for multivalent counterions.
- **Physical Review E Manuscript (`doc/paper_pre.tex` & `doc/paper_pre.pdf`)**:
  - Corrected minor typographical mole fractions in Table III for System 5C ($x_{-1} = 0.09, x_{-2} = 0.02$).
  - Overhauled Section IV.D to present the complete multi-system 2PT analysis with Table IV and updated Figure 6.
  - Recompiled manuscript to `doc/paper_pre.pdf` (12 pages, clean compilation).

---

## [V1.4] - 2026-10-21

### Added
- **Species-Resolved Partial Molar Excess Entropies & Partial Molar Volumes in HNC**:
  - Implemented numerical temperature- and density-derivatives of excess chemical potentials $\mu_i^{\text{ex}}$ across all four multicomponent mixtures.
  - Computed the canonical constant-volume partial excess entropy $s_{i,V}^{\text{ex}} = -\left(\frac{\partial \mu_i^{\text{ex}}}{\partial T}\right)_V$, the constant-pressure partial molar excess entropy $\bar{S}_{i,P}^{\text{ex}} = -\left(\frac{\partial \mu_i^{\text{ex}}}{\partial T}\right)_P$, and partial molar volumes $\bar{V}_i$.
  - Demonstrated that weakening secondary solvent dispersion ($\varepsilon_2 = 1.25 \to 0.75$) yields a massive surge in its partial molar excess entropy ($\Delta \bar{S}_{2,P}^{\text{ex}} = +0.465\,k_B$, rising from $-1.580\,k_B$ to $-1.115\,k_B$).
  - Quantified electrostriction and valence ordering: partial molar volumes decrease systematically with ionic charge ($\bar{V}_{\text{solv}} \approx 1.446 > \bar{V}_{|z|=1} \approx 1.439 > \bar{V}_{|z|=2} \approx 1.418\,\sigma^3$), while partial molar excess entropies become increasingly negative ($\bar{S}_{\text{solv}}^{\text{ex}} = -1.357\,k_B > \bar{S}_{|z|=1}^{\text{ex}} = -1.423\,k_B > \bar{S}_{|z|=2}^{\text{ex}} = -1.619\,k_B$).
  - Updated Physical Review E manuscript (`doc/paper_pre.tex`) with new Section II.E (2PT theoretical formulation), Section IV.C (Partial molar excess entropies & electrostriction, Table III), Section IV.D (2PT analysis & benchmark against TI, Table IV, Figure 6), and Conclusions; recompiled `doc/paper_pre.pdf` (11 pages).
- **Two-Phase Thermodynamics (2PT) Molecular Dynamics Framework (`lammps_study/two_pt_analysis/`)**:
  - Prepared and executed an NVT simulation with high-frequency velocity dumping on an NVIDIA RTX PRO 4500 GPU.
  - Implemented `two_pt_calculator.py` parsing trajectories, evaluating species velocity autocorrelation functions (VACF), vibrational densities of states (DoS), fluidicity parameter $f$, and 2PT entropy components.
  - Benchmarked 2PT against MD Thermodynamic Integration (TI): Green-Kubo diffusion quantitatively matches Einstein MSD diffusion ($D = 0.2152$ vs $0.2148$), and 2PT excess entropy accurately captures the constant-volume partial excess entropy ($s_{i,V}^{\text{ex}} = -5.369\,k_B$ vs HNC $-5.268\,k_B$, within $1.8\%$).

---

## [V1.3] - 2026-10-21

### Added
- **Species-Resolved Partial Pair Entropy ($s_i^{(2)}$) Analysis**:
  - Investigated the origin of the decoupling between divalent ion diffusion and bulk excess entropy.
  - Implemented the Green-Wallace species-resolved partial pair entropy:
    $$s_i^{(2)} = - 2\pi \rho \sum_{j=1}^p x_j \int_0^\infty \left[ g_{ij}(r) \ln g_{ij}(r) - (g_{ij}(r)-1) \right] r^2 dr$$
    evaluated directly from MD radial distribution functions $g_{ij}(r)$ across all multicomponent fluids.
  - Demonstrated that while bulk excess entropy is dominated by the majority solvent matrix ($89$--$94\%$), multivalent ions ($|z|=2$) decouple into deeply negative structural caging ($s_i^{(2)} \approx -1.72$ to $-1.76\,k_B$).
  - Identified contact ion-pairing ($g_{+-}(r_{\text{contact}})$ rising from $4.93$ in 3C to $5.25$ in 4C $\varepsilon_2=0.75$) as the microscopic origin of suppressed divalent ion mobility in weakened solvent mixtures.
  - Established a universal master curve for all neutral solvent species obeying $D_i^* \propto \exp(B s_i^{(2)})$ with $B = +1.05$ ($R^2 = 0.91$).

### Changed
- **Figure 5 Redesign (`doc/figures/rosenfeld_scaling.png`)**:
  - Re-engineered Figure 5 into a two-panel publication layout:
    - **Panel (a)**: Net fluid Rosenfeld scaling for the controlled solvent-variation series (System 3C, 4C $\varepsilon_2=1.25$, and 4C $\varepsilon_2=0.75$, sharing identical $I=0.12$ and $6\%$ divalent ions), yielding a strictly positive scaling exponent $B = +1.05$ (MD TI) and $B = +1.40$ (HNC), alongside System 5C as an independent compositional baseline ($I=0.19$, $80\%$ solvent).
    - **Panel (b)**: Species-resolved scaling of reduced diffusion coefficients $D_i^*$ against partial pair entropy $s_i^{(2)}$.
  - Resolved regression slope discrepancy ($B = -0.14$ caused by mixing disparate compositional families in a single unconstrained fit).
- **Physical Review E Manuscript (`doc/paper_pre.tex` & `doc/paper_pre.pdf`)**:
  - Updated Section IV.D and Conclusions with detailed physical discussion on species-resolved pair entropy scaling, contact ion pairing, and multivalent counterion transport decoupling.

---

## [V1.2] - 2026-10-21

### Added
- **4-Component Fluid with Weakened Solvent Dispersion ($\varepsilon_2 = 0.75$)**:
  - Investigated chemical heterogeneity and dispersion asymmetry by introducing a less attractive solvent ($89\%$ Solvent 1 with $\varepsilon_1=1.0$, $5\%$ Solvent 2 with $\varepsilon_2=0.75$, $3\%$ Cation 2+, $3\%$ Anion 2-, $I=0.12$).
  - Performed complete HNC density sweep down to $\rho^* = 0.01$ and confirmed that reducing solvent dispersion weakens local caging and increases excess entropy ($s_{\text{ex}} = -2.0551\,k_B$ at $\rho^* = 0.70$, compared to $-2.0635\,k_B$ for 3C and $-2.0703\,k_B$ for $\varepsilon_2=1.25$).
  - Configured and executed a complete batch of 15 GPU-accelerated LAMMPS MD simulations across the supercritical isotherm ($T^*=2.00, \beta^*=0.50$) on NVIDIA RTX PRO 4500 Blackwell GPUs (`ladon28`, `ladon29`).
- **Thermodynamic Integration & Entropy Across Four Multicomponent Fluids**:
  - Generalized MD thermodynamic integration pipeline (`thermodynamic_integration_md.py`) to process 3C (Option 1, $I=0.12$), 4C ($\varepsilon_2=1.25$), 4C ($\varepsilon_2=0.75$), and 5C ($I=0.19$).
  - Demonstrated quantitative concordance between MD and HNC across all 4 systems, validating the entropic ordering: $s_{\text{ex}}(4\text{C}, \varepsilon_2=1.25) < s_{\text{ex}}(3\text{C}) < s_{\text{ex}}(4\text{C}, \varepsilon_2=0.75)$.
- **Self-Diffusion Dynamics & Rosenfeld Scaling**:
  - Evaluated high-precision self-diffusion coefficients from 150,000-step production trajectories at $\rho^* = 0.70$ ($R^2 > 0.999$).
  - Uncovered direct solvent-dispersion coupling: Solvent 2 diffuses slower ($D = 0.2032$) when more attractive ($\varepsilon_2 = 1.25$) and significantly faster ($D = 0.2348$) when less attractive ($\varepsilon_2 = 0.75$), driving an overall increase in net fluid mobility ($D_{\text{net}} = 0.2194$).
  - Confirmed quantitative validity of Rosenfeld's excess entropy scaling relation across all four fluids with positive scaling exponent $B \approx +0.75$.
- **Physical Review E Manuscript Update**:
  - Expanded `doc/paper_pre.tex` and recompiled `doc/paper_pre.pdf` incorporating the four-system comparative study.
  - Updated all 5 multi-panel publication figures in `doc/figures/` (`compressibility_and_ti.png`, `entropy_comparison.png`, `delta_entropy.png`, `msd_comparison.png`, `rosenfeld_scaling.png`).

---

## [V1.1] - 2026-10-21

### Added
- **Molecular Dynamics Simulation Suite (`lammps_study/`)**:
  - GPU-accelerated LAMMPS simulation framework for multicomponent charged Lennard-Jones mixtures ($N = 2000$ particles) along the supercritical isotherm ($T^* = 2.00, \Gamma = 1.0$) across 15 density points ($\rho^* \in [0.01, 0.70]$).
  - SLURM batch submission pipelines targeting NVIDIA RTX PRO 4500 Blackwell GPUs on cluster nodes `ladon28` and `ladon29`.
  - Full long-range Coulomb electrostatics via PPPM with adaptive real-space cutoff ($r_c = 8.0\,\sigma$ for $\rho^* \le 0.05$ to collapse reciprocal mesh overhead, and $r_c = 2.5\,\sigma$ for $\rho^* \ge 0.10$).
  - Full trajectory recording at $\rho^* = 0.70$ in both NetCDF (`traj.nc`) and custom LAMMPS (`traj.lammpstrj`) formats.
- **Thermodynamic Integration from Molecular Dynamics**:
  - Time-averaged equations of state ($P, Z, U_{\text{ex}}$) integrated via cubic splines along the isotherm.
  - Evaluation of excess Helmholtz free energy $A_{\text{ex}} / (N k_B T)$ and total excess entropy $s_{\text{ex}} / (N k_B)$ from simulation.
  - Excellent agreement with HNC internal energy ($|\Delta u_{\text{ex}}| < 0.2\%$ across $\rho^* \in [0.10, 0.35]$) and thermodynamic offsets ($\Delta s_{\text{ex}}^{\text{MD}} = -0.0223\,k_B$ vs $\Delta s_{\text{ex}}^{\text{HNC}} = -0.0192\,k_B$).
- **Diffusion Dynamics & Rosenfeld Excess Entropy Scaling**:
  - Center-of-mass corrected Mean Squared Displacements (MSD) over 500,000 production steps ($1000\,\tau$).
  - High-precision self-diffusion coefficients ($R^2 > 0.999$) for net fluid, total divalent ions ($|z|=2$), cations, and anions.
  - Quantitative demonstration of Rosenfeld's scaling hypothesis: enhanced structural ordering at higher ionic strength retards diffusion, yielding strictly positive scaling exponents $B_{\text{net}} \approx +0.69$--$+1.10$ and $B_{\text{an}} \approx +4.67$.
- **Physical Review E Manuscript (`doc/paper_pre.tex` & `doc/paper_pre.pdf`)**:
  - Complete RevTeX 4.2 publication manuscript formatted for Physical Review E.
  - Integrated theory (multicomponent OZ, Allnatt screening, Morita-Hiroike excess free energy functional, exact Euler relation $Z_{\text{diff}} \equiv Z_{\text{virial}}$, Green-Wallace pair entropy $s^{(2)}$).
  - Multi-panel publication figures in `doc/figures/`:
    - `compressibility_and_ti.png`: Equation of state $Z(\rho^*)$ and free energy $A_{\text{ex}}(\rho^*)$ comparing MD against HNC virial and free-energy routes.
    - `entropy_comparison.png`: Excess entropy $s_{\text{ex}}(\rho^*)$ and pair entropy $s^{(2)}(\rho^*)$ along the isotherm.
    - `delta_entropy.png`: Entropic offset $\Delta s_{\text{ex}}(\rho^*)$ and $\Delta s^{(2)}(\rho^*)$ comparing MD and HNC.
    - `msd_comparison.png`: Center-of-mass corrected MSD (ballistic crossover and linear Einstein diffusive regime).
    - `rosenfeld_scaling.png`: Rosenfeld excess entropy scaling plots for net fluid and divalent ions.

---

## [V1.0] - 2026-10-20

### Added
- **Licensing**:
  - Adopted GNU General Public License v3.0 (GPLv3), aligned with `~/trj_analysis`.
- **Modern Modular Architecture**:
  - Ported legacy monolithic `HNCR-nCLJ.f90` to modern Fortran 2008/2018.
  - Split into 10 cohesive modules under `src/`:
    - `types_mod.f90`: Standard floating-point precisions and physical constants.
    - `params_mod.f90`: Runtime problem dimensions and parameters.
    - `grid_mod.f90`: Dynamic radial ($r$) and reciprocal ($q$) grid infrastructure.
    - `fft_mod.f90`: Fast Sine/Cosine transforms with precomputed bit-reversal lookup tables.
    - `arrays_mod.f90`: Dynamic heap-allocated arrays for correlations and Jacobians.
    - `linalg_mod.f90`: BLAS (`DGEMV`) and LAPACK (`DGETRF`, `DGETRI`) linear algebra wrappers.
    - `potentials_mod.f90`: Cut-and-shifted Lennard-Jones potential, virials, and Debye screening.
    - `initial_guess_mod.f90`: Analytical MSA for charged spheres and Percus--Yevick HS initializations.
    - `thermo_mod.f90`: Internal energy, Helmholtz free energy, chemical potentials, and structure factors.
    - `hnc_solver_mod.f90`: Coupled Broyles mixing and Newton-Raphson nonlinear solver with Kinoshita-Harada continuation.
    - `main.f90`: Command-line interface and I/O coordination.
- **Build Infrastructure & Optimization**:
  - `Makefile` supporting GNU (`gfortran`), Intel (`ifx`/`ifort`), and NVIDIA (`nvfortran`).
  - Automated library detection (OpenBLAS, FlexiBLAS, Intel MKL).
  - Validation test target (`make test`).
  - Delivers a 25% CPU performance boost over legacy implementation.
- **Entropy Comparison Study (`EntropyComparison/`)**:
  - Comparative investigation of excess entropy $s_{\text{ex}}$ and two-body pair entropy $s^{(2)}$ between:
    - **System 3C**: Solvent ($x_1=0.80$) + pure symmetric 2:2 electrolyte ($q = [0, +2, -2]$, $x = [0.80, 0.10, 0.10]$, ionic strength $I = 0.40$).
    - **System 5C**: Solvent ($x_1=0.80$) + mixed 1:1 and 2:2 electrolytes ($q = [0, +1, +2, -1, -2]$, $x = [0.80, 0.05, 0.04, 0.09, 0.02]$, ionic strength $I = 0.19$).
  - Systems perfectly matched in packing fraction ($\eta \approx 0.3665$ at $\rho^* = 0.70$), solvent volume fraction, total ion volume fraction, and Lennard-Jones parameters ($\sigma = 1.0, \epsilon = 1.0, r_{\text{cut}} = 2.5$) at electrostatic coupling $\gamma = 1.0$.
  - Selected supercritical isotherm at $T^* = 2.00$ ($\beta^* = 0.50 > \beta_c^* \approx 0.76$) ensuring $\chi_T^{-1} > 0$ throughout, safely avoiding any vapor-liquid spinodal down to $\rho^* = 0.01$.
  - Automated pipeline script `EntropyComparison/scripts/run_study.sh`.
  - Python post-processing and thermodynamic integration suite `EntropyComparison/scripts/thermodynamic_integration.py`:
    - Evaluates $A_{\text{ex}} / N k_B T = \int_0^\rho \frac{Z-1}{\rho'} d\rho'$ via cubic splines.
    - Demonstrates that the 3C system (higher ionic strength) exhibits greater structural ordering ($\Delta s_{\text{ex}} \approx +0.019 \, k_B$ for 5C relative to 3C), with $\sim 58\%$ of this entropy difference captured by the two-body pair term $\Delta s^{(2)}$.
    - Generates publication plots (`entropy_comparison.png`, `delta_entropy.png`, `compressibility_and_ti.png`) and tabular data (`comparison_data.csv`, `comparison_summary.txt`).
- **Two-Body Structural Pair Entropy ($s^{(2)}$)**:
  - Implemented exact multicomponent Green--Wallace pair entropy integral:
    $$\frac{s^{(2)}}{N k_B} = - \frac{1}{2} \rho \sum_{i,j} x_i x_j \int \left[ g_{ij}(r) \ln g_{ij}(r) - (g_{ij}(r) - 1) \right] d\mathbf{r}$$
  - Evaluated via Simpson quadrature across radial grid points with analytical core handling ($g \to 0 \implies \text{integrand} \to 1$).
  - Added `s2/nk` as an aligned column in `Thermohncr.dat`.
- **Excess & Total Entropy Outputs**:
  - Excess entropy: $s_{\text{ex}} / N k_B = (u_{\text{ex}} - a_{\text{ex}}) / k_B T$.
  - Ideal-gas mixture entropy: $s_{\text{id}} / N k_B = \frac{5}{2} - \ln \rho^* - \sum_i x_i \ln x_i$.
  - Total classical entropy: $s_{\text{tot}} = s_{\text{id}} + s_{\text{ex}}$.
- **LaTeX Theoretical Documentation** (`doc/entropy_theory.tex` & `doc/entropy_theory.pdf`):
  - In-depth technical note explaining the occurrence of negative classical continuous entropy at low temperatures / high ionic strengths.
  - Connection to continuous Gibbs entropy divergence ($\lim_{T \to 0} S = -\infty$).
  - Quantification of the omitted quantum thermal de Broglie wavelength offset in reduced units ($-3\ln(\Lambda/\sigma) \approx +8 \text{ to } +12 \, k_B$).
  - Mathematical properties of $s^{(2)}$ (negative definiteness, core behavior, asymptotic decay) and relation to Rosenfeld scaling.
- **5-Component Test Suite (`5C/`)**:
  - Sample dataset `5C/datosLJ.hncr` with 5 charged species ($q = [0, +1, +2, -1, -2]$, $x = [0.80, 0.05, 0.04, 0.09, 0.02]$) across temperatures $\beta^* = 0.10 \dots 1.10$.
  - Verified thermodynamic consistency between virial compressibility $Z_{\text{virial}}$ and free-energy/chemical-potential compressibility $Z_a$ ($< 0.5\%$ deviation).

### Fixed
- **Thermodynamic Consistency in Free Energy ($a_{\text{ex}}^{\text{HNC}}$)**:
  - Corrected core radius reference in [`src/thermo_mod.f90`](file:///home/e.lomba/EntropicFLuids/src/thermo_mod.f90): in adapting from the hard-sphere RPM solver, the core cavity volume ($a_{\text{SR}}$), chemical potential hard-core shift, and electrostatic surface terms were still evaluated using the hard-sphere diameter $\sigma(j)$ rather than the soft-core cutoff boundary $R_{\text{core}}(j) = R(ncore(j)) \approx 0.6 \, \sigma(j)$.
  - Updating these terms removes an artificial $+1.642 \, \rho$ linear slope discrepancy in $a_{\text{ex}}^{\text{HNC}}(\rho)$, restoring exact internal thermodynamic consistency:
    - 5-point numerical differentiation $Z = 1 + \rho \frac{\partial a_{\text{ex}}}{\partial \rho}$ now reproduces the virial compressibility factor $Z_{\text{virial}}$ with **$0.01\% - 0.05\%$ accuracy** across all densities.
    - Thermodynamic integration along the isotherm reproduces the direct HNC free energy to within **$< 0.01 \, k_B T$**.

### Changed
- **Arbitrary Multicomponent Scaling ($NSP \ge 2$)**:
  - Generalized from previous 4-species limit to arbitrary $NSP$ (tested up to 5 components).
  - Analytical closed-form electroneutrality solver for the last two species ($x_{NSP-1}$ and $x_{NSP}$) from $NSP-2$ user-specified mole fractions.
  - Dynamically weighted Bhatia-Thornton structure factors ($S_{cc}$) and effective solvent-ion direct correlation functions ($c_{\text{ion-ion}}$, $c_{\text{solvent-ion}}$).
  - Expanded internal linear algebra workspace in `linalg_mod.f90` from $6\times 6$ to $32\times 32$.
- **Unified Output Formatting**:
  - Formatted `Thermohncr.dat` with dynamic headers aligning to exact character widths across arbitrary $NSP$.
  - Converted `RDFhncr.dat`, `skhnc.dat`, and `lsol.hncr` formatting to standard Fortran 2008 unlimited repeat specifiers `*(...)` to eliminate line wrapping when $NIT > 10$.
- **Build System**:
  - Auto-detection of OpenBLAS paths and compiler runtime rpaths in `Makefile`.
  - Added cleaning of stale `.mod` files in project root.
