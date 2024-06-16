#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=NPT_solv
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=./logs/NVT_nonp/%x-%j.out
#SBATCH --error=./logs/NVT_nonp/%x-%j.err

srun python3 alchemical_osoa_NPT.py -m 'cyclopentanol'
srun python3 alchemical_osoa_NVT_non_periodic.py -m 'cyclopentanol'