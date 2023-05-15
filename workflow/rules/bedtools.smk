#Version 2.30.0
#Python version: 3.8.12

####################################################################
# Intersect GIAB HC variants vcf with Sequencing Platform's BED file
####################################################################
rule giab_intersect_bed:
    input:
        GIAB_HC_VARIANTS_FILE
    output:
        GIAB_INTERSECT_FILE
    conda:
        f'{ENVS_DIR}/bedtools.yaml'
    threads: 2 
    params:
        bedfile=get_bed
    priority: 0 #did not use before
    shell:
        f"bedtools intersect -a {{input}} -b {{params.bedfile}} -u -header > {{output}}"

###########################################################
# Intersect VarDict vcf with Sequencing Platform's BED file
###########################################################
rule sample_intersect_bed:
    input:
        VARDICT_FILE
    output:
        SRA_INTERSECT_FILE
    conda:
        f'{ENVS_DIR}/bedtools.yaml'
    threads: 2  
    params:
        bedfile=get_bed
    priority: 0 #did not use before
    shell:
        f'bedtools intersect -a {{input}} -b {{params.bedfile}} -u -header > {{output}}'

