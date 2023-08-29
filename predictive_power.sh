#!/bin/bash
#
#SBATCH --job-name=predictive_power
#SBATCH --output=/home/klam/seq_artifacts/results/logs/predictive_power9_joined_LR_feats_bsF.log
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00


# Purge and load modules:
#module purge
#module load sratoolkit-3.0.2
#module load fastqc
#module load bwa
#module load samtools
#module load picard
#module load intel-python3

eval "$(conda shell.bash hook)"
conda activate jupyterlab

python3 predictive_power.py