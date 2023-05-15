#Version: 1.10.1
#Python version: 3.8.12

#####################################################
# Joined quality control of all raw and trimmed reads
#####################################################
rule multiqc_quality:
    input:
        ALL_FASTQC_FILES,
        ALL_TRIMMED_FASTQC_FILES
    conda:
        f'{ENVS_DIR}/multiqc.yaml'
    output:
        MULTIQC_FILE
    shell:
        f'multiqc -f -n {MULTIQC_NAME} -c {MULTIQC_CONFIG_FILE} -o {MULTIQC_DIR} {{input}}'

        