# outputs for this rule will be stored in the resources directory
# all other rule outputs will be in the results/outputs directory

#Version: sratoolkit-3.0.2 (HPC module)

# prefetch and fasterq-dump used together is the fastest way to extract FASTQ-files from SRA-accessions
###################################
# Downloads all necessary sra files
###################################
rule prefetch:
    output:
        SAMPLE_SRA_FILE
#    conda:
#        f'{ENVS_DIR}/sratoolkit.yaml'
    threads: 4
    priority: 11
    shell:
        f'prefetch --option-file {SRA_LIST_TXT} -O {SRA_DIR}'

##############################################
# Extracts reads from the downloaded sra files
##############################################
rule fasterq_dump:
    input:
        SAMPLE_SRA_FILE
    output:
        SAMPLE_F_FILE,
        SAMPLE_R_FILE
#    conda:
#        f'{ENVS_DIR}/sratoolkit.yaml'
    threads: 8
    priority: 10 
    shell:
        f'''
        fasterq-dump {{input}} --split-files --outdir {DATA_DIR}
        gzip {DATA_DIR}/*.fastq
        '''

