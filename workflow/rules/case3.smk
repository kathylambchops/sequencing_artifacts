#Version 2.30.0
#Python version: 3.8.12

####################################################
# Bedtools coverage to distinguish case 3.1 from 3.2
####################################################
rule bedtools_coverage:
    input:
        a=CASE3_REGIONS_BED,
        b=COORD_SORTED_BAM_FILE
    output:
        CASE3_COVERAGE_BED
    conda:
        f'{ENVS_DIR}/bedtools.yaml'
    threads: 2
    priority: 0
    shell:
        f'bedtools coverage -a {{input.a}} -b {{input.b}} -d > {{output}}'


#a BAM/BED/GFF/VCF file “A”. Each feature in A is compared to B in search of overlaps. Use “stdin” if passing A with a UNIX pipe.
#b One or more BAM/BED/GFF/VCF file(s) “B”. Use “stdin” if passing B with a UNIX pipe. NEW!!!: -b may be followed with multiple databases and/or wildcard (*) character(s).
#d Report the depth at each position in each A feature. Positions reported are one based. Each position and depth follow the complete A feature.



# Version:2.26.10 (HPC tool)

#####################################################
# Add read group to BAM file for GATK haplotypecaller
#####################################################
# This tool enables the user to assign all the reads in the #INPUT to a single new read-group.
rule picard_add_readgroup:
    input:
        COORD_SORTED_BAM_FILE
    output:
        RG_COORD_SORTED_BAM_FILE
#    conda:
#        f'{ENVS_DIR}/picard.yaml'
    priority: 0 
    threads: 2
    params:
        platform=get_platform
    shell:
        f'''
        java -Xmx5G -jar /opt/picard/picard.jar AddOrReplaceReadGroups \
        I={{input}} \
        O={{output}} \
        RGLB=lib1 \
        RGPL={{params.platform}} \
        RGPU=unit1 \
        RGSM=23 \
        CREATE_INDEX=true
        '''

#I Input file (BAM or SAM or a GA4GH url).
#O Output file (BAM or SAM).
#RGLB Read-Group library (required - arbitrary)
#RGPL Read-Group platform (e.g. ILLUMINA, SOLID)
#RGPU Read-Group platform unit (eg. run barcode) (required - arbitrary)
#RGSM Read-Group sample name (required - arbitrary)
#CREATE_INDEX here because this BAM is the one we will use for variant calling - GATK haplotypecaller

 

#Version 4.2.4.1 (HPC tool)

###########################################################
# GATK haplotypecaller to confirm case 3.2 call matches REF 
###########################################################
rule GATK_haplotypecaller:
    input:
        l=CASE3_REGIONS_BED,
        i=RG_COORD_SORTED_BAM_FILE
    output:
        GATK_FILE
    conda:
         f'{ENVS_DIR}/gatk.yaml'
    threads: 2
    priority: 0
    shell:
        f'''
        gatk --java-options "-Xmx5G" HaplotypeCaller \
        -R {GENOME_FILE} \
        -I {{input.i}} \
        -O {{output}} \
        -ERC GVCF \
        -L {{input.l}} \
        -ip 100
        '''
        
#I BAM/SAM/CRAM file containing reads
#O File to which variants should be written
#R Reference sequence file (needs .dict and .fai file in same directory)
#L One or more genomic intervals over which to operate. Run the variant calling with a targeted regions 
#ip Amount of padding (in bp) to add to each interval you are including. Use this to add padding to the intervals specified using -L. 
#   This is typically used to add padding around targets when analyzing exomes.
#   For example, '-L 1:100' with a padding value of 20 would turn into '-L 1:80-120'.
#   make sure to give sufficient padding around the targeted sites (100 bp on each side)
