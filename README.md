# Identifying and Characterizing Sequencing Artifacts in Next-Generation Sequencing Data using Machine Learning Methods

Next Generation Sequencing (NGS) coupled with liquid biopsies is a non-invasive way to detect low cancer mutations, which have low allelic frequencies. This technique can be easily repeated without any risk or side effects, making it an appealing approach.  The issue at hand is that during NGS, artifactual variants arise from DNA library preparation methods and errors in the NGS platforms. These artifacts can be mistaken as true variants and affect the accuracy of variant calling techniques. Because of this, it is vital to be able to distinguish between real variants and sequencing artifacts. This project builds a bioinformatics pipeline to process Whole Exome Sequencing (WES) data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) database. Comparison of variant calls to Genome In a Bottle (GIAB) high confidence (HC) regions allows differentiation between artifacts and non-artifacts. These labels are used to train and test supervised machine learning (ML) models such as Random Forest, Extra Trees, and XGBoost. The left flanking region (LSEQ) and right flanking region (RSEQ) of each single nucleotide polymorphism (SNP) in the variant calling file were extracted. Metrics such as the content of each base, each kmer of size 4 and its count, the largest homopolymer size, the largest palindrome size, and the largest hairpin loop size were computed to be used as features in model building.

## Pair-end DNA-sequencing pipeline
This pipeline was built to run on San Jose State University (SJSU) High Performance Computer (HPC).
To run the pipeline using slurm:
```
sbatch run_snakemake.sh
```

## Setup
Snakemake and mamba must be installed.  To create a conda environment with the required packages:
```
conda create -n bfx_env python=3.11
conda activate bfx_env
conda install -c conda-forge mamba
mamba install snakemake -c conda-forge -c bioconda
```
The specific versions used in this project were:
- Python (v._3.11.0_)
- Snakemake (v._7.21.0_)
- Mamba (v._1.3.0_)

## Configs
Update the following configuration files:
- `config/`
  - **[config.yaml](https://github.com/kathylambchops/sequencing_artifacts/blob/main/config/config.yaml)**: Contains settings for the Snakemake workflow (_used in [`workflow/Snakefile`](https://github.com/kathylambchops/sequencing_artifacts/blob/main/workflow/Snakefile)_)
  - **[multiqc_config.yaml](https://github.com/kathylambchops/sequencing_artifacts/blob/main/config/multiqc_config.yaml)**: Contains settings for the MultiQC report (_used in [`config/config.yaml`](https://github.com/kathylambchops/sequencing_artifacts/blob/main/config/config.yaml)_)
- `slurm/`
  - **[config.yaml](https://github.com/kathylambchops/sequencing_artifacts/blob/main/slurm/config.yaml)**: Contains settings for running the pipeline using Slurm (_used in [`run_snakemake.sh`](https://github.com/kathylambchops/sequencing_artifacts/blob/main/run_snakemake.sh)_)
  - **[cluster.yaml](https://github.com/kathylambchops/sequencing_artifacts/blob/main/slurm/cluster.yaml)**: Contains settings for job submissions (_used in [`slurm/config.yaml`](https://github.com/kathylambchops/sequencing_artifacts/blob/main/slurm/config.yaml)_)


## Results
The outputs folder, specified by **outputs_dir** in [`config/config.yaml`](https://github.com/kathylambchops/sequencing_artifacts/blob/main/config/config.yaml), should have the following directory structure:
```
/path/to/outputs
├── bedtools
├── bwa
├── fastqc
├── logs
├── multiqc
├── picard
├── vardict
└── trimmomatic
```

## Resources
All resources, including data and reference genome, must be placed inside the **resources** folder. Filepaths are written relative to it when specified in [`config/config.yaml`](https://github.com/kathylambchops/sequencing_artifacts/blob/main/config/config.yaml).
The directory structure is as follows:
```
/path/to/resources
├── adapters
├── bedfiles
├── data
├── GIAB
├── hg38
└── sra
```


