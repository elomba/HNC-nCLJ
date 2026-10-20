# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
