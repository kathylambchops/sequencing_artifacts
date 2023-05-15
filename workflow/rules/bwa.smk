#Version: 0.7.15-r1140 (HPC tool)

##############################################
# Map adapter trimmed reads to referene genome
##############################################
rule bwa_map:
    input:
        f=TRIMMED_F_PAIRED_FILE,
        r=TRIMMED_R_PAIRED_FILE
    output:
        MAPPED_BAM_FILE
    priority: 9 
#    conda:
#        f'{ENVS_DIR}/bwa.yaml'
    threads: 10
    shell:
        f"time bwa mem -t {{threads}} {GENOME_FILE} {{input.f}} {{input.r}} | samtools view -Sb -@ {{threads}} > {{output}}"
# do - > if doesnt work





