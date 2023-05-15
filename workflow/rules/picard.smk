# Version:2.26.10 (HPC tool)

#############################
# Sort aligned reads by query
#############################
rule picard_sort_query:
    input:
        MAPPED_BAM_FILE
    output:
        QUERY_SORTED_BAM_FILE
#    conda:
#        f'{ENVS_DIR}/picard.yaml'
    priority: 8 
    shell:
        f'java -Xmx5G -jar /opt/picard/picard.jar SortSam -I {{input}} -O {{output}} -SORT_ORDER queryname -TMP_DIR {TMP_DIR}'

#I The SAM, BAM or CRAM file to sort.  Required.
#O  The sorted SAM, BAM or CRAM output file.   Required
#SORT_ORDER Sort order of output file


#######################################################
# Remove duplicates from the query sorted aligned reads
#######################################################
rule picard_dedupe:
    input:
        QUERY_SORTED_BAM_FILE
    output:
        b=DEDUPED_BAM_FILE,
        l=DEDUPED_LOG_FILE
#    conda:
#        f'{ENVS_DIR}/picard.yaml' 
    priority: 8 
    shell:
        f'java -Xmx5G -jar /opt/picard/picard.jar MarkDuplicates -I {{input}} -O {{output.b}} -M {{output.l}} -TMP_DIR {TMP_DIR} \
        -REMOVE_DUPLICATES true -ASSUME_SORT_ORDER queryname'

#I(string) One or more input SAM or BAM files to analyze. Must be coordinate sorted. Required
#O(file) The output file to write marked records to Required.
#M(file) File to write duplication metrics to  Required.
#CREATE_INDEX Whether to create an index when writing VCF or coordinate sorted BAM output.
#ASSUME_SORT_ORDER assume that the input file has this order


###############################################
# Sort aligned deduplicated reads by coordinate
###############################################
rule picard_sort_coord:
    input:
        DEDUPED_BAM_FILE
    output:
        COORD_SORTED_BAM_FILE
#    conda:
#        f'{ENVS_DIR}/picard.yaml'
    priority: 7 
    shell:
        f'java -Xmx5G -jar /opt/picard/picard.jar SortSam -I {{input}} -O {{output}} -SORT_ORDER {SORT_ORDER} -TMP_DIR {TMP_DIR} -CREATE_INDEX true'
        
#CREATE_INDEX here because this sorted BAM is the one we will use for variant calling - Vardict