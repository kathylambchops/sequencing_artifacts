#!/bin/bash
#
#SBATCH --job-name=lazypredict
#SBATCH --output=/home/klam/seq_artifacts/results/logs/lazypredict2.log
#
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=20000


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

python3 lazypredict_models.py