#!/bin/bash
#SBATCH --job-name=s2_md_pilot
#SBATCH --partition=ladon_v4_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/sif/users/e.lomba/EntropicFLuids/MDPairEntropyAssessment/logs/pilot_%j.out
#SBATCH --error=/sif/users/e.lomba/EntropicFLuids/MDPairEntropyAssessment/logs/pilot_%j.err
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /sif/users/e.lomba/EntropicFLuids
python3 MDPairEntropyAssessment/test_pair_entropy.py
python3 MDPairEntropyAssessment/extract.py --index 0 --limit 20
python3 MDPairEntropyAssessment/extract.py --index 21 --limit 20
