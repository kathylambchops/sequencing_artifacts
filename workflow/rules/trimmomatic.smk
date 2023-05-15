#Version: 0.39 (HPC tool)

##############################
# Adapter readthrough trimming
##############################
rule trimmomatic_trim:
    input:
        f=SAMPLE_F_FILE,
        r=SAMPLE_R_FILE
    output:
        fp=TRIMMED_F_PAIRED_FILE,
        fup=TRIMMED_F_UNPAIRED_FILE,
        rp=TRIMMED_R_PAIRED_FILE,
        rup=TRIMMED_R_UNPAIRED_FILE
#    conda:
#        f'{ENVS_DIR}/trimmomatic.yaml'
    priority: 10 
    threads: 8
    shell:
        f'time java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE {{input.f}} {{input.r}} \
        {{output.fp}} {{output.fup}} {{output.rp}} {{output.rup}} SLIDINGWINDOW:4:20 ILLUMINACLIP:{ADAPTERS_FILE}:2:30:10'
        
