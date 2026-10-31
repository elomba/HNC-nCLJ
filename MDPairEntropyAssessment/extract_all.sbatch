#!/bin/bash
#SBATCH --job-name=s2_md_extract
#SBATCH --partition=ladon_v4_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --array=0-23%12
#SBATCH --output=/sif/users/e.lomba/EntropicFLuids/MDPairEntropyAssessment/logs/extract_%A_%a.out
#SBATCH --error=/sif/users/e.lomba/EntropicFLuids/MDPairEntropyAssessment/logs/extract_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /sif/users/e.lomba/EntropicFLuids
python3 MDPairEntropyAssessment/extract.py --index "$SLURM_ARRAY_TASK_ID"
python3 MDPairEntropyAssessment/analyze.py --index "$SLURM_ARRAY_TASK_ID"
python3 MDPairEntropyAssessment/extract.py --index "$SLURM_ARRAY_TASK_ID" --ideal
python3 MDPairEntropyAssessment/analyze.py --index "$SLURM_ARRAY_TASK_ID" --ideal
