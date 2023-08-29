# Identifying and Characterizing Sequencing Artifacts in Next-Generation Sequencing Data using Machine Learning Methods

Next Generation Sequencing (NGS) coupled with liquid biopsies is a non-invasive way to detect low cancer mutations, which have low allelic frequencies. This technique can be easily repeated without any risk or side effects, making it an appealing approach.  The issue at hand is that during NGS, artifactual variants arise from DNA library preparation methods and errors in the NGS platforms. These artifacts can be mistaken as true variants and affect the accuracy of variant calling techniques. Because of this, it is vital to be able to distinguish between real variants and sequencing artifacts. This project builds a bioinformatics pipeline to process Whole Exome Sequencing (WES) data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) database. Comparison of variant calls to Genome In a Bottle (GIAB) high confidence (HC) regions allows differentiation between artifacts and non-artifacts. These labels are used to train and test supervised machine learning (ML) models such as Random Forest, Extra Trees, and XGBoost. The left flanking region (LSEQ) and right flanking region (RSEQ) of each single nucleotide polymorphism (SNP) in the variant calling file were extracted. Metrics such as the content of each base, each kmer of size 4 and its count, the largest homopolymer size, the largest palindrome size, and the largest hairpin loop size were computed to be used as features in model building.

The full paper can be found here: [Characterizing Sequencing Artifacts](https://scholarworks.sjsu.edu/cgi/viewcontent.cgi?article=2279&context=etd_projects).

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
/path/to/results
├── logs
├── outputs
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


## Machine Learning pipeline
The `final.ipynb` notebook was developed to import the `get_metrics.py` module built for this project. All of the SRA VCF files were converted into individual Pandas data frames with the LSEQs and RSEQs. The flanking region metrics were computed and populated into the data frames. Some Artifact data exploration was also done here. The python module `get_metrics.py` was created to store all the functions that determines artifacts and computes flanking region metrics.  More information about how each function works can be found here: [Characterizing Sequencing Artifacts](https://scholarworks.sjsu.edu/cgi/viewcontent.cgi?article=2279&context=etd_projects).

A final data frame(s) used for ML is created using `ML_df_creation.ipynb`

The lazypredict package was used to build a lot of basic models without any parameter tuning to identify which top models would be used for additional ML tasks. You should have a final data frame you plan on using the lazypredict script for. Edit `lazy_predict.py` with the correct file that contains your final data frame, then run the script.
To run the script using slurm:
```
sbatch lazypredict.sh
```

Predictive power and feature importance were determined using `predictive_power.py`.  This trains 1000 random forest models. Run as many times as needed to get enough coverage, in my case I ran 10 instances to get a total of 10000 models trained before determining feature importance.  Edit `predictive_power.py` with the correct file that contains your final data frame and change the output file name for each instance you run this script.
To run the script using slurm:
```
sbatch predictive_power.sh
```

Join the predictive power subfiles together (if you ran multiple instances) so you are left with a single file before investigating predictive power. `LazypredictResults_InvestigatePredictivePower.ipynb` looks at the output from running `lazypredict.py` and/or `lazy_predict.sh`.  It also investigates feature importance.

Finally, `ModelEvaluation.ipynb` evaluates model performance of all features, just the LSEQ/RSEQ features I calculate from `get_metrics.py`, only LSEQ/RSEQ features with 60% or greater predictive power based on `predictive_power.py`, and from LSEQ/RSEQ feature clustering based on a package called VarClusHi. 
