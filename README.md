# HNCR-nCLJ: Hypernetted Chain Solver for Multicomponent Charged Lennard-Jones Fluids

**HNCR-nCLJ** is a high-performance Fortran 2008/2018 code that solves the multicomponent Ornstein-Zernike (OZ) integral equation with the Hypernetted Chain (HNC) closure for mixtures of charged Lennard-Jones particles (soft-core ions) for an arbitrary number of species (e.g. 2, 3, 4, 5+ components).

The code implements **Allnatt's renormalized OZ formulation** to handle long-range Coulombic forces without divergence, coupled with the **Labík–Malijevský–Vonka** Newton-Raphson iterative procedure and the **Kinoshita–Harada** continuation strategy across temperature and density sweeps.

---

## Features

- **Soft-Core Ionic Model**: Cut-and-shifted Lennard-Jones (12-6) core plus Coulombic electrostatic interactions with arbitrary per-pair parameters ($\sigma_{ij}$, $\epsilon_{ij}$, $r_{\text{cut},ij}$) and ionic valencies ($q_i$).
- **Arbitrary Multicomponent Support**: Handles 2 to 5+ species with automatic, closed-form electroneutrality and stoichiometric mole fraction determination.
- **Renormalized Integral Equations**: Allnatt's formulation subtracting Yukawa-screened Coulomb tails in both real and Fourier space.
- **Fast & Robust Convergence**:
  - Inner **Newton-Raphson** loop on the direct correlation function over short-range grid points.
  - Outer **Broyles** mixing loop with an adaptive mixing parameter $\alpha$.
  - **Kinoshita–Harada continuation**: The inverted Jacobian matrix is preserved across steps, drastically speeding up thermodynamic state sweeps.
- **Optimized Numerics**:
  - Precomputed trigonometric and bit-reversal lookup tables in Fast Sine/Cosine Transforms (FFTSIN/FFTCOS).
  - High-performance BLAS (`DGEMV`) and LAPACK (`DGETRF`, `DGETRI`) matrix algebra with auto-detection of OpenBLAS, FlexiBLAS, and Intel MKL.
- **Thermodynamic Consistency & Structure**:
  - Internal energy (Coulomb and LJ contributions).
  - Helmholtz free energy $A / N k_B T$ via $k$-space integration.
  - **Excess entropy** $S_{\text{ex}} / N k_B = (U_{\text{ex}} - A_{\text{ex}}) / N k_B T$.
  - **Total entropy** $S / N k_B = S_{\text{id}} / N k_B + S_{\text{ex}} / N k_B$ including ideal-gas mixture mixing entropy.
  - **Two-body pair entropy** $S^{(2)} / N k_B = -\frac{1}{2}\rho \sum_{i,j} x_i x_j \int \left[ g_{ij}(r)\ln g_{ij}(r) - (g_{ij}(r)-1) \right] d\mathbf{r}$ evaluated via Simpson quadrature.
  - Virial and free-energy compressibility factors ($Z_v$ and $Z_a$).
  - Reduced inverse isothermal compressibility $\chi_T^{-1}$.
  - Excess chemical potentials $\mu_i / k_B T$.
  - Partial pair correlation functions $g_{ij}(r)$ and structure factors $S_{ij}(q)$.
  - Bhatia-Thornton structure factors ($S_{NN}(q)$, $S_{cc}(q)$) and concentration fluctuation eigenvalues ($\lambda_1, \lambda_2$).

---

## Theory and Model

### 1. Pair Interaction Potential
The interaction potential between particles of species $i$ and $j$ separated by distance $r$ is:
$$u_{ij}(r) = u_{ij}^{\text{LJ}}(r) + u_{ij}^{\text{Coulomb}}(r)$$

where the truncated and shifted Lennard-Jones potential is:
$$u_{ij}^{\text{LJ}}(r) = \begin{cases} 4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right] - u_{ij}^{\text{cut}}, & r \le r_{\text{cut}, ij} \\ 0, & r > r_{\text{cut}, ij} \end{cases}$$
$$u_{ij}^{\text{cut}} = 4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{\text{cut}, ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{\text{cut}, ij}}\right)^6\right]$$

and the Coulomb potential is:
$$u_{ij}^{\text{Coulomb}}(r) = \frac{q_i q_j e^2}{\epsilon r}$$

### 2. Allnatt Renormalization
To eliminate the $1/r$ divergence in the Fourier transforms of ionic fluids, the potential and direct correlation functions are decomposed into short-range ($S$) and long-range ($L$) parts:
$$c_{ij}(r) = c_{ij}^{S}(r) + c_{ij}^{L}(r), \quad s_{ij}(r) = s_{ij}^{S}(r) + s_{ij}^{L}(r)$$
The long-range electrostatic tails are screened using a Yukawa-like term with inverse Debye length $\kappa_D = \sqrt{\frac{4\pi}{\beta \epsilon} \sum_k \rho_k q_k^2}$.

### 3. Closure Relation
The Hypernetted Chain (HNC) closure in renormalized form reads:
$$c_{ij}(r) = \exp\left[-u_{ij}^{\text{LJ}}(r) + h_{ij}(r) - c_{ij}(r)\right] - 1 - (h_{ij}(r) - c_{ij}(r))$$

### 4. Entropy Formulation and Two-Body Pair Contribution
The code evaluates three complementary entropy quantities:
1. **Excess Entropy**:
   $$\frac{s_{\text{ex}}}{N k_B} = \frac{u_{\text{ex}} - a_{\text{ex}}}{N k_B T} < 0$$
   quantifying the configurational ordering induced by excluded-volume repulsion and electrostatic ion pairing.
2. **Two-Body Structural Entropy (Green--Wallace Expansion)**:
   $$\frac{s^{(2)}}{N k_B} = - 2\pi \rho \sum_{i=1}^{NSP} \sum_{j \le i} (2 - \delta_{ij}) x_i x_j \int_0^\infty \left[ g_{ij}(r) \ln g_{ij}(r) - (g_{ij}(r) - 1) \right] r^2 dr \le 0$$
   which is strictly negative definite and captures $85\text{--}95\%$ of the total excess entropy.
3. **Total Classical Entropy**:
   $$\frac{s}{N k_B} = \left[ \frac{5}{2} - \ln \rho^* - \sum_{i=1}^{NSP} x_i \ln x_i \right] + \frac{s_{\text{ex}}}{N k_B}$$
   *Note on negative total entropy*: In classical continuous statistical mechanics ($\lim_{T\to 0} S = -\infty$) and reduced units ($\Lambda = \sigma = 1$), the omitted quantum thermal de Broglie wavelength contribution $-3\ln(\Lambda/\sigma) \approx +8 \text{ to } +12 \, k_B$ shifts the baseline downward. At low temperatures ($\beta^* \ge 0.9$), $|s_{\text{ex}}| > s_{\text{id}}$, causing reduced total entropy to become negative while remaining entirely physical. For full derivations, see [`doc/entropy_theory.pdf`](file:///home/e.lomba/EntropicFLuids/doc/entropy_theory.pdf) or [`doc/entropy_theory.tex`](file:///home/e.lomba/EntropicFLuids/doc/entropy_theory.tex).

---

## Project Structure

```
.
├── Makefile                # Unified build system (GNU, Intel, NVIDIA)
├── README.md               # Documentation and usage guide
├── CHANGELOG.md            # Detailed version history
├── datosLJ.hncr            # Default input parameters (3-component charged mixture)
├── 5C/                     # Dedicated 5-component test suite
│   ├── datosLJ.hncr        # 5-species input dataset
│   └── Thermohncr.dat      # Converged 5-component thermodynamic output
├── EntropyComparison/      # 3C vs. 5C excess entropy comparative study & TI
│   ├── system_3C/          # 3-component (solvent + 2:2 salt) inputs
│   ├── system_5C/          # 5-component (solvent + 1:1 & 2:2 salt) inputs
│   ├── scripts/            # Pipeline (run_study.sh) & TI script (thermodynamic_integration.py)
│   ├── results/            # Density sweep datasets, comparison summary, and plots
│   └── README.md           # Methodology, physical parameters, and detailed analysis
├── doc/                    # Theoretical notes and LaTeX documentation
│   ├── entropy_theory.tex  # Comprehensive entropy & reduced units technical note
│   └── entropy_theory.pdf  # Compiled PDF document
├── examples/               # Example input datasets
│   ├── datosLJ_3comp.hncr
│   └── datosLJ_2comp_nocharge.hncr
├── baseline_results/       # Reference outputs for numerical validation
└── src/
    ├── types_mod.f90       # Standard precision and physical constants
    ├── params_mod.f90      # Global parameters and species limits
    ├── grid_mod.f90        # Radial (r) and reciprocal (q) grid management
    ├── fft_mod.f90         # Fast Sine/Cosine Transforms with bit-reversal tables
    ├── arrays_mod.f90      # Dynamic allocation of correlation functions and matrices
    ├── linalg_mod.f90      # LAPACK/BLAS wrappers (invmat, matvec) & symmetric solvers
    ├── potentials_mod.f90  # Truncated LJ potential, virial, and screening functions
    ├── initial_guess_mod.f90# Analytical MSA and Percus-Yevick starting estimates
    ├── thermo_mod.f90      # Energy, pressure, chemical potentials, Bhatia-Thornton S(q), S_ex, S^(2)
    ├── hnc_solver_mod.f90  # Coupled Broyles / Newton-Raphson iterative engine
    └── main.f90            # Main driver and CLI entry point
```

---

## Prerequisites and Compilation

### Prerequisites
- Modern Fortran compiler (`gfortran` $\ge 9$, `ifx`, or `ifort`).
- BLAS and LAPACK libraries (e.g., OpenBLAS, FlexiBLAS, or Intel MKL).

### Building with `make`

To build the optimized release binary (`bin/hncr_nclj`):
```bash
make
```

To build with full debugging checks (bounds checking, tracebacks, floating-point traps):
```bash
make MODE=debug
```

To specify a custom compiler or linear algebra library:
```bash
# Using Intel oneAPI (ifx + MKL)
make FC=ifx LIBS="-qmkl"

# Using FlexiBLAS
make FC=gfortran LIBS="-lflexiblas"
```

To run the automated validation test:
```bash
make test
```

To clean build artifacts:
```bash
make clean      # Removes obj/ directory
make distclean  # Removes obj/, bin/, and calculation outputs
```

---

## Running the Code

Execute the compiled binary:
```bash
./bin/hncr_nclj [input_file]
```
If no input file is passed as an argument, the program defaults to `datosLJ.hncr` in the current working directory.

---

## Input File Specification (`datosLJ.hncr`)

| Card # | Variable(s) | Description |
|:---|:---|:---|
| **1** | `title` | Title / description string (up to 80 chars) |
| **2** | `NR, DR, tdiff, pmixt, tmax, amx` | Grid points ($2^{\text{NU}}$), grid step $\Delta r$, convergence tolerance, Broyles mixing parameter, max CPU seconds, mixing error threshold |
| **3** | `NUM, diffnr` | Number of grid points for Newton-Raphson scheme ($\le ND$), NR tolerance |
| **4** | `NSP` | Number of species ($NSP \ge 2$) |
| **5** | `ssig(1:NSP, 1:NSP)` | Collision diameters $\sigma_{ij}$ matrix ($NSP \times NSP$) |
| **6** | `eeps(1:NSP, 1:NSP)` | Lennard-Jones well depths $\epsilon_{ij}$ matrix ($NSP \times NSP$) |
| **7** | `qq(1:NSP)` | Ionic charges $q_i$ in elementary units |
| **8** | `nrho, ntr` | Number of total densities and temperatures |
| **9** | `densid(1:nrho)` | Total number densities $\rho^* = (N/V)\sigma^3$ |
| **10** | `tk(1:ntr)` | Reduced temperatures $T^* = k_B T / \epsilon$ |
| **11** | `compsi` | Logical flag (`.true.` = compute MSA initial guess, `.false.` = read `solin.hncr`) |
| **12** | `XXc(1:NSP-2)` | Independent mole fractions (species $1 \dots NSP-2$; remaining 2 determined by electroneutrality) |
| **13** | `rct(1:NSP, 1:NSP)` | Interaction cutoff radii $r_{\text{cut}, ij}$ matrix ($NSP \times NSP$) |
| **14** | `gama` | Solvent dielectric / electrostatic coupling scaling parameter |

---

## Output Files

| File | Content |
|:---|:---|
| `Salida.hncr` | Detailed log of iterations, residuals, and timing per thermodynamic state |
| `Thermohncr.dat` | Thermodynamic summary table: density $\rho$, temperature $\beta^*$, pressure $Z$, energies ($U_{\text{tot}}, U_{\text{coulomb}}$), free energy $A/Nk_BT$, excess entropy $S_{\text{ex}}/Nk_B$, total entropy $S/Nk_B$, two-body entropy $S^{(2)}/Nk_B$, compressibility $\chi_T^{-1}$, chemical potentials $\mu_i$, and Bhatia-Thornton eigenvalues $\lambda_1, \lambda_2$ |
| `RDFhncr.dat` | Radial distribution functions $g_{ij}(r)$ for all pairs as a function of $r$ |
| `skhnc.dat` | Structure factors $S_{ij}(q)$, total number $S_{NN}(q)$, and concentration-concentration $S_{cc}(q)$ |
| `lsol.hncr` | Converged short-range series function $r \cdot s_{ij}(r)$ (can be used to restart with `compsi=.false.`) |

---

## References & Historical Lineage

1. **A.R. Allnatt**, *Renormalization of the Ornstein-Zernike Equation for Ionic Solutions*, Mol. Phys. **8**, 533 (1964).
2. **S. Labík, A. Malijevský, and P. Vonka**, *An efficient iterative method for solving integral equations of fluids*, Mol. Phys. **56**, 709 (1985).
3. **M. Kinoshita and M. Harada**, *An improved numerical method for solving the hypernetted chain equation for fluid mixtures*, Mol. Phys. **65**, 599 (1988).
4. **E. Lomba**, *On the solution of the reference hypernetted chain equation for simple fluids*, Mol. Phys. **68**, 87 (1989).
5. **E. Lomba and J.S. Hoye**, Trondheim (1991).
6. **E. Lomba and C. Martín**, Madrid (2007).
7. **E. Lomba**, Ferrol (2020).

---

## License
Academic / Open Source. See included files for original notices.
