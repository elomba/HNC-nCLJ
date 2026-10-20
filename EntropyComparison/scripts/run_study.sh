#!/usr/bin/env bash
#===============================================================================
# Script: run_study.sh
# Description: Automated execution of the 3-Component vs. 5-Component
#              excess entropy comparison and density continuation sweeps.
#===============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${BASE_DIR}/.." && pwd)"
BIN="${ROOT_DIR}/bin/hncr_nclj"

echo "========================================================================"
echo " Starting Entropy Comparison Study: 3-Component vs 5-Component Mixtures"
echo "========================================================================"
echo "Executable: ${BIN}"

if [[ ! -x "${BIN}" ]]; then
    echo "Executable not found. Building project..."
    (cd "${ROOT_DIR}" && make)
fi

mkdir -p "${BASE_DIR}/results"

#-------------------------------------------------------------------------------
# 1. System 3C: Solvent + z=2 Electrolyte (q = [0, +2, -2])
#-------------------------------------------------------------------------------
echo ""
echo ">>> [1/4] System 3C: Initializing at rho=0.70 up to beta*=0.50..."
cd "${BASE_DIR}/system_3C"
"${BIN}" input_init_rho07.hncr
cp -f lsol.hncr solin.hncr
cp -f Thermohncr.dat "${BASE_DIR}/results/Thermo_3C_init_rho07.dat"

echo ""
echo ">>> [2/4] System 3C: Density continuation sweep (rho = 0.70 -> 0.01)..."
"${BIN}" input_density_sweep.hncr
cp -f Thermohncr.dat "${BASE_DIR}/results/Thermo_3C_density_sweep.dat"

#-------------------------------------------------------------------------------
# 2. System 5C: Solvent + 1:1 and 2:2 Electrolytes (q = [0, +1, +2, -1, -2])
#-------------------------------------------------------------------------------
echo ""
echo ">>> [3/4] System 5C: Initializing at rho=0.70 up to beta*=0.50..."
cd "${BASE_DIR}/system_5C"
"${BIN}" input_init_rho07.hncr
cp -f lsol.hncr solin.hncr
cp -f Thermohncr.dat "${BASE_DIR}/results/Thermo_5C_init_rho07.dat"

echo ""
echo ">>> [4/4] System 5C: Density continuation sweep (rho = 0.70 -> 0.01)..."
"${BIN}" input_density_sweep.hncr
cp -f Thermohncr.dat "${BASE_DIR}/results/Thermo_5C_density_sweep.dat"

#-------------------------------------------------------------------------------
# 3. Post-Processing & Thermodynamic Integration
#-------------------------------------------------------------------------------
echo ""
echo ">>> Running thermodynamic integration and comparative analysis..."
cd "${BASE_DIR}/scripts"
python3 thermodynamic_integration.py

echo ""
echo "========================================================================"
echo " Entropy Comparison Study Completed Successfully!"
echo " Results stored in: ${BASE_DIR}/results/"
echo "========================================================================"
