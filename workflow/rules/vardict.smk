#Version: 1.8.3=hdfd78af_0

##############################
# Quality control of raw reads
##############################
rule variant_calling:    
    input:
        COORD_SORTED_BAM_FILE 
    output:
        VARDICT_FILE
    conda:
        f'{ENVS_DIR}/vardict.yaml'
    threads: 14
    priority: 6 
    shell:
        f'''
        start=$(date +%s)
        vardict-java \
        -th 32 \
        -G {GENOME_FILE} \
        -f {ALLELE_FREQ} \
        -b {{input}} \
        -c 1 -S 2 -E 3 -g 4 {GIAB_BED_FILE} \
        | teststrandbias.R | var2vcf_valid.pl \
        -E -f {ALLELE_FREQ} > {{output}}
        end=$(date +%s)
        runtime_s=$(echo $(( end - start )))
        sec_per_min=60
        sec_per_hr=3600
        runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
        echo "total run time(m): $runtime_m"
        runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
        echo "total run time(h): $runtime_h"
        '''
                
                