#Version: 0.11.9
#Python version: 3.8.12

##############################
# Quality control of raw reads
##############################
rule fastqc_quality:
    input:
        SAMPLE_F_FILE,
        SAMPLE_R_FILE
    output:
        FASTQC_F_FILE,
        FASTQC_R_FILE
    conda:
        f'{ENVS_DIR}/fastqc.yaml'
    priority: 5 
    threads: 2
    shell:
        f'fastqc --noextract -t {{threads}} -o {FASTQC_DIR} {{input}}'

##########################################
# Quality control of adapter trimmed reads
##########################################
rule fastqc_quality_trimmed:
    input:
        TRIMMED_F_PAIRED_FILE,
        TRIMMED_R_PAIRED_FILE
    output:
        TRIMMED_FASTQC_F_FILE,
        TRIMMED_FASTQC_R_FILE
    conda:
        f'{ENVS_DIR}/fastqc.yaml'
    priority: 5 
    threads: 2
    shell:
        f'fastqc --noextract -t {{threads}} -o {FASTQC_DIR} {{input}}'