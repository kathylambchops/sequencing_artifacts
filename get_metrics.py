# module to determine artifacts and get metrics

import pandas as pd
import numpy as np

def vcf_to_df(vcf_file, extract_flank_seqs=False): 
    """ Converts vcf file into a pandas dataframe. Extracts CHROM, POS, REF, ALT. Set extract_flank_seqs=True
        if vcf is generated from vardict and LSEQ, RSEQ are desired. AF, HIAF, HICNT, VD, SN, ADJAF, VARBIAS 
        are also extracted because these features have high percent accuracy for model training.
    """
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith('#CHROM'):
                header_names = line.strip().split()
                header_names[0] = header_names[0][1:]
                break
        read_lines = vcf.readlines()
        # list containing each row as 1 giant string

    
    # split the giant string into columns -- now we have a 2D list
    create_columns = [row.strip().split() for row in read_lines]

    final = []
    
    # extracts LSEQ and RSEQ for SRR vcf file if needed, in addition to CHROM, POS, REF, ALT
    # also extract TYPE (SNV, Insertion, Deletion, Complex)
    if extract_flank_seqs == True:
        for row in create_columns:
            # CHROM, POS, REF, ALT, LSEQ, RSEQ, AF, HIAF, HICNT, VD, SN, ADJAF, VARBIAS
            final.append([row[0], row[1], row[3], row[4], 
                          row[7].split(';')[24][5:], 
                          row[7].split(';')[25][5:], 
                          row[7].split(';')[4][3:], 
                          row[7].split(';')[16][5:],
                          row[7].split(';')[22][6:],
                          row[7].split(';')[3][3:],
                          row[7].split(';')[15][3:],
                          row[7].split(';')[17][6:],
                          row[7].split(';')[7][8:]]) 
        df = pd.DataFrame(data=final, columns=header_names[0:2]+header_names[3:5]+
                          [row[7].split(';')[24][:4], 
                          row[7].split(';')[25][:4], 
                          row[7].split(';')[4][:2], 
                          row[7].split(';')[16][:4],
                          row[7].split(';')[22][:5],
                          row[7].split(';')[3][:2],
                          row[7].split(';')[15][:2],
                          row[7].split(';')[17][:5],
                          row[7].split(';')[7][:7]])
                          #['LSEQ','RSEQ', 'AF', 'HIAF', 'HICNT', 'VD', 'SN', 'ADJAF', 'VARBIAS']
        return df
    
    # For general use vcfs - extracts CHROM, POS, REF, ALT    
    for row in create_columns:
        final.append([row[0], row[1], row[3], row[4]]) # CHROM, POS, REF, ALT                     
    df = pd.DataFrame(data=final, columns=header_names[0:2]+header_names[3:5])   
    return df

    
def filter_for_snvs(df):
    """ Function filters for SNVs only by ensuring the string length of 'REF' and 'ALT' is 1.
    """    
    return df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]

def snp_freq_not1(df):
    """ Takes df and finds positions where ALT has multiple variant calls (ex: C, CAA). If any of the calls for that position are 
        single nucleotides, then that position is considered to contain a SNP. Returns df of these types of SNPs
    """
    alt_series = df.ALT.loc[lambda alt_col: [',' in alt for alt in alt_col]]
    final_consensus = []
    for j in alt_series:
        split_alt = j.split(',')
        row_consensus = False
        for sub_alt in split_alt:
            if len(sub_alt)==1:
                row_consensus = True
                break
        final_consensus.append(row_consensus)
    
    df2_idx = alt_series.loc[final_consensus].index.tolist()    
    df2 = df.loc[df2_idx]
    df2 = df2[df2['REF'].str.len() == 1]    
    return df2

def join_datasets(srr_df, giab_df):
    """ Merge SRR and GIAB pandas df on 'CHROM' and 'POS'
    """
    joined = pd.merge(
        # Merge the SRR and GIAB dataframes...
        srr_df,
        giab_df,

        # Using their indexes (which are CHROM and POS)...
        left_index=True,
        right_index=True,

        # Using suffixes according to the dataset names...
        suffixes=("_SRR", "_GIAB"),

        # Taking all rows from both sides (even if they don't have a match)
        how='outer'
    )
    return joined


def extract_artifacts_column(joined_dfs):
    """ Extract Artifacts and differentiate between 4 different cases
    """
    is_same_alt = joined_dfs.ALT_SRR == joined_dfs.ALT_GIAB # Variant is same and at the same location 
    in_srr = joined_dfs.ALT_SRR.notnull() # Includes all the variants, except NaN's
    in_giab = joined_dfs.ALT_GIAB.notnull()
    
    # Series with the mutiple variant calls -- exclude NaN values
    multi = joined_dfs[(joined_dfs.ALT_GIAB.str.len() != 1) & in_giab & in_srr] 
    # If one of the ALT_GIAB multiple calls is the same as the ALT_SRR call then it is not an artifact
    same = [(alt_srr in alt_giab.split(',')) for alt_srr, alt_giab in zip(multi.ALT_SRR, multi.ALT_GIAB)]
    # Update is_same_alt boolean values with the special multi cases
    is_same_alt.loc[multi.index] = same    

    non_artifacts = in_srr & in_giab & is_same_alt # Case 1 -- true variants
    different_variant = in_giab & in_srr & ~is_same_alt # Case 2 -- ~ bitwise not
    not_in_srr = in_giab & ~in_srr # Case 3
    not_in_giab = in_srr & ~in_giab # Case 4

    print("Case 1:", non_artifacts.sum())
    print("Case 2:", different_variant.sum())
    print("Case 3:", not_in_srr.sum())
    print("Case 4:", not_in_giab.sum())
    
    is_artifact = different_variant | not_in_giab # | is or -- considers case 2 and 4 only
    # ignore case 3 as it is difficult to determine 3.1 from 3.2 right now
    #is_artifact = different_variant | not_in_giab | not_in_srr
    
    return is_artifact


def content(dna_list, flank):
    """Takes a list of sequences and the string L or R for which flanking sequence.
        Returns a list of dictionaries storing the content of A, T, C, G
    """
    content_list = []
    for seq in dna_list:
        a = seq.upper().count('A') / len(seq) * 100
        t = seq.upper().count('T') / len(seq) * 100
        g = seq.upper().count('G') / len(seq) * 100
        c = seq.upper().count('C') / len(seq) * 100
        a,t,g,c = round(a,2), round(t,2), round(g,2), round(c,2) 
        content_list.append({flank+'_A':a, flank+'_T':t, flank+'_G':g, flank+'_C':c})
    return content_list


def find_kmer(dna_list, flank, kmer_size=4):
    """Takes a list of DNA sequences and an integer for kmer_size. 
        Finds all kmers and their counts.  Returns a list of dictionaries
        Each dict contains all kmers and their counts for each seq
    """ 
    # kmer_dict_list holds a dict for each lseq/rseq
    # kmer_dict_list is -- [{seq1_4mer},{seq2_4mer},{seq3_4mer}...etc]
   
    kmer_dict_list = [] 
    for seq in dna_list: #dna_list=['GGGGGTGGTGG']
        kmer_list = []
        # kmer lists for each seq
        # Any sequence of length L will contain (L - k + 1) total k-mers
        for i in range(len(seq)-kmer_size+1):
            kmer = seq[i:i+kmer_size]
            kmer_list.append(kmer)
        # kmer_count() returns a dictionary
        kmer_dict = kmer_count(kmer_list, flank)
        kmer_dict_list.append(kmer_dict)            
    
    return kmer_dict_list

# helper function for find_kmer()
def kmer_count(kmer_list, flank):
    """Given a list of kmers, find the count of each unique one. Return dictionary where
        kmer is key and count is value.
    """
    kmer_dict = {}
    for kmer in kmer_list:   #kmer_list ['GGGG', 'GGGG', 'GGGT', 'GGTG', 'GTGG', 'TGGT', 'GGTG', 'GTGG']
        kmer = flank+'_'+kmer
        if kmer not in kmer_dict:
            kmer_dict[kmer] = 1    
        else:
            kmer_dict[kmer] += 1           
    return kmer_dict


# homopolymer 3 or more of a base (yes/no) -- do for each base
def homopolymer(dna_list, flank, min_polymer_size=3):
    """Takes a list of DNA sequences and an integer for minimum polymer size. Calls find_kmer
       function and determines if kmers are homopolymers. Does this for all kmer sizes from min_polymer_size
       to length of the sequence.  Returns list of dictionaries showing largest homopolymer size for each base"""
    
    # call find_kmer function which was previously defined to find kmer of size min_polymer_size and their counts (counts irrelevant here)
    kmer_dict_list = find_kmer(dna_list, flank, min_polymer_size)
    # returns something like [{'R_GGG': 4, 'R_GGT': 2, 'R_GTG': 2, 'R_TGG': 1}, ...]
    
    homopoly_list = []
    for d in kmer_dict_list:
        homo_poly_d = {}
        # if kmer contains all the same bases then it is a homopolymer of size k
        if flank+'_'+('A'*min_polymer_size) in d.keys():
            homo_poly_d[flank+'_HOMO_POLY_A'] = min_polymer_size
        else:
            homo_poly_d[flank+'_HOMO_POLY_A'] = 0
        if flank+'_'+('T'*min_polymer_size) in d.keys():
            homo_poly_d[flank+'_HOMO_POLY_T'] = min_polymer_size
        else:
            homo_poly_d[flank+'_HOMO_POLY_T'] = 0
        if flank+'_'+('G'*min_polymer_size) in d.keys():
            homo_poly_d[flank+'_HOMO_POLY_G'] = min_polymer_size
        else:
            homo_poly_d[flank+'_HOMO_POLY_G'] = 0
        if flank+'_'+('C'*min_polymer_size) in d.keys():
            homo_poly_d[flank+'_HOMO_POLY_C'] = min_polymer_size
        else:
            homo_poly_d[flank+'_HOMO_POLY_C'] = 0
        homopoly_list.append(homo_poly_d)

    # test for homopolymers of up to the length of the LSEQ/RSEQ 
    for polymer_size in range(min_polymer_size + 1, len(dna_list[0]) + 1):  # 4, 21
        kmer_dict_list = find_kmer(dna_list, flank, polymer_size)
        #print(kmer_dict_list)
        
        # j is used to index homopolymer_list so each dictionary in that list can be modified 
        j = 0

        for d in kmer_dict_list:

            homo_poly_d2 = {}
            
            if flank+'_'+('A'*polymer_size) in d.keys():
                homo_poly_d2[flank+'_HOMO_POLY_A'] = polymer_size
            # no larger homopolymer found, then keep the previous value that is already stored in the initial dictionary
            else:
                homo_poly_d2[flank+'_HOMO_POLY_A'] = homopoly_list[j][flank+'_HOMO_POLY_A']
            if flank+'_'+('T'*polymer_size) in d.keys():
                homo_poly_d2[flank+'_HOMO_POLY_T'] = polymer_size
            else:
                homo_poly_d2[flank+'_HOMO_POLY_T'] = homopoly_list[j][flank+'_HOMO_POLY_T']
            if flank+'_'+('G'*polymer_size) in d.keys():
                homo_poly_d2[flank+'_HOMO_POLY_G'] = polymer_size
            else:
                homo_poly_d2[flank+'_HOMO_POLY_G'] = homopoly_list[j][flank+'_HOMO_POLY_G']
            if flank+'_'+('C'*polymer_size) in d.keys():
                homo_poly_d2[flank+'_HOMO_POLY_C'] = polymer_size
            else:
                homo_poly_d2[flank+'_HOMO_POLY_C'] = homopoly_list[j][flank+'_HOMO_POLY_C']
            
            # homo_poly_d2 contains updated data for each polymer size -- update the previous dictionary with new dictionary
            homopoly_list[j].update(homo_poly_d2)
            
            # increment j to access next dictionary in homo_poly_list
            j += 1
            
    return homopoly_list   


# palindrome - form complimentary with itself? -- sliding window to see largest palindrome size.

# Original single strand                5’ - ACCTAGGT - 3’
# Complement                            3’ - TGGATCCA - 5’
# Reverse complement == OG strand       5’ - ACCTAGGT - 3’   PALINDROME

# Original single strand                5’ - TTCTCTGC - 3’
# Complement                            3’ - AAGAGACG - 5’
# Reverse complement != OG strand       5’ - GCAGAGAA - 3’  NOT A PALINDROME

def palindrome(dna_list, flank):
    """Takes a list of DNA sequences and flank. Calls find_kmer function and determines if kmers are palindromes. 
       Does this for all kmer sizes from 3 to length of the sequence. Find the reverse complement of each 
       kmer sequence by calling complement() function. If reverse complement is identical to original 
       kmer sequence, it is a palindrome. Returns list containing largest palindrome size for each seq.     
    """
    # ordered array filled with 0's representing palindrome size for each seq -- if palindrome exists, update 0 value with largest palindrome size
    palindrome_list = [0] * len(dna_list)
    
    # test for palindromes of size 4 up to the length of the LSEQ/RSEQ 
    for kmer_size in range(4, len(dna_list[0]) + 1):  # 4, 21
        kmer_dict_list = find_kmer(dna_list, flank, kmer_size)
        # returns something like [{'R_GGGA': 4, 'R_GGTA': 2, 'R_GTGA': 2, 'R_TGGA': 1}, ...]
        
        # j used to index palindrome_list
        j = 0
        
        for d in kmer_dict_list:
            # paldrome_d has kmer seq as key and size
            #palindrome_d = {}
            # if kmer is a palindrome then the sequence has palindrome of size kmer_size
            for kmer in d:
                kmer = kmer.lstrip(flank+'_') # gives kmer without the L_ and R_ tags
                comp = complement(kmer)
                rev_complement = comp[::-1]
                if rev_complement == kmer and palindrome_list[j] < kmer_size:
                    # update size of palindrome for that seq based on index
                    palindrome_list[j] = kmer_size 
            # update j so we can move on to next seq
            j += 1
        
    return palindrome_list


# helper function for palindrome()
def complement(seq):
    """Takes a DNA sequence. Returns the complement of the sequence"""
    complement = seq.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()
    return(complement)


# hairpin loop -- does it fold back on itself to form a paired double helix
# optimal loop length is 4-8 bases
# length of LSEQ and RSEQ is 20 bp
# min stem size is 4

def hairpin(dna_list, loop_size = 4):
    """Takes a list of DNA sequences and an integer for loop_size.  Optimal loop_size is between 4 and 8.
       Determines if each sequence contains a hairpin loop or not. If flanking regions of the loop are 
       reverse complement of eachother, then it folds back on itself to form a hairpin loop. This portion
       is called the stem.  
    """
    # hairpin loop: ---stem1---loop---stem2---
   
    # all possible stem lengths
    stem_lengths = [4,5,6,7,8]
    hairpin_list = []  # list of boolean values shows if each sequence has a hairpin loop or not
    
    for seq in dna_list:
        # contains the sequence of the loop if a hairpin is formed
        loop_list = []
        for stem_length in stem_lengths:      
            count = 0
            # Any sequence of length L will contain (L - k + 1) total k-mers
            # treat the stem1---loop---stem2 unit as the 'kmer' 
            kmer_size = stem_length*2 + loop_size
            
            for i in range(len(seq)-kmer_size+1):  # (L - k + 1)
                kmer = seq[i:i+kmer_size]
                stem1 = kmer[:stem_length]
                loop = kmer[stem_length:stem_length+loop_size]
                stem2 = kmer[len(kmer)-stem_length:]
                count += 1
                #print(stem1,loop,stem2)
                # if reverse complement of stem1 == stem2 then it is a haripin because it will fold back on itself
                if complement(stem1)[::-1] == stem2:
                    loop_list.append(loop)
                    
            #print(count)
            #print(loop_list) 
           
        if len(loop_list) > 0:
            hairpin_list.append(True)
        else:
            hairpin_list.append(False)            
        
    return hairpin_list


def main(SRR_path):
    """ Takes path to SRR VCF file and calls on all previous functions to return a complete df containing all metrics
    """
    # columns in df: ['CHROM', 'POS', 'REF', 'ALT', 'LSEQ', 'RSEQ', 'AF', 'HIAF', 'HICNT', 'VD', 'SN', 'ADJAF', 'VARBIAS']
    df = vcf_to_df(SRR_path, extract_flank_seqs=True)
    
    # create a separate df for lseq and rseq. This way all lseq metric columns will follow lseq
    # and all rseq metric columns will follow rseq
    # the final df will merge the two together -- this is easier than using insert() over and over
    
    df_lseq = df[['CHROM','POS','REF','ALT','AF','HIAF','HICNT','VD','SN','ADJAF','VARBIAS','LSEQ']]
    df_rseq = df[['CHROM','POS','REF','ALT','RSEQ']]

    lseq = df['LSEQ']  # list of left flanking sequences
    rseq = df['RSEQ']  # list of right flanking sequences

    ##############################
    # call content() function with list of lseqs's to find content of each bp
    lseq_content = content(lseq, 'L')

    # call content() function with list of rseqs's to find content of each bp
    rseq_content = content(rseq, 'R')

    # add lseq_content to lseq df and rseq_content to rseq df
    df_lseq = pd.concat([df_lseq, pd.DataFrame(lseq_content)], axis=1)
    df_rseq = pd.concat([df_rseq, pd.DataFrame(rseq_content)], axis=1)

    ##############################
    # call find_kmer() function with list of lseqs's to find all 4mers 
    lseq_kmers = find_kmer(lseq, 'L', 4)

    # call find_kmer() function with list of rseqs's to find all 4mers 
    rseq_kmers = find_kmer(rseq, 'R', 4)

    # add lseq_kmer to lseq df and rseq_kmer to rseq df
    df_lseq = pd.concat([df_lseq, pd.DataFrame(lseq_kmers).fillna(0).astype(int)], axis=1)
    df_rseq = pd.concat([df_rseq, pd.DataFrame(rseq_kmers).fillna(0).astype(int)], axis=1)

    ##############################
    # call homopolymer() function with list of lseqs's to determine if there are homopolymers of size 3 
    lseq_homo_poly = homopolymer(lseq, 'L', 3)

    # call homopolymer() function with list of rseqs's to determine if there are homopolymers of size 3  
    rseq_homo_poly = homopolymer(rseq, 'R', 3)

    # add lseq_homo_poly to lseq df and rseq_homo_poly to rseq df
    df_lseq = pd.concat([df_lseq, pd.DataFrame(lseq_homo_poly)], axis=1)
    df_rseq = pd.concat([df_rseq, pd.DataFrame(rseq_homo_poly)], axis=1)

    ##############################
    # call palindrome() function with list of lseqs's to determine if the seq's are palindromes or not
    lseq_palindrome = palindrome(lseq, 'L')

    # call palindrome() function with list of rseqs's to determine if the seq's are palindromes or not  
    rseq_palindrome = palindrome(rseq, 'R')

    # add lseq_palindrome to lseq df and rseq_palindrome to rseq df
    df_lseq = pd.concat([df_lseq, pd.DataFrame(lseq_palindrome, columns=['L_PALINDROME'])], axis=1)
    df_rseq = pd.concat([df_rseq, pd.DataFrame(rseq_palindrome, columns=['R_PALINDROME'])], axis=1)

    ##############################
    # call hairpin() function with list of lseqs's to determine if there is a hairpin loop of size 4-8
    lseq_hairpin_matrix = [hairpin(lseq, i) for i in range(4,9)] 

    # call hairpin() function with list of rseqs's to determine if there is a hairpin loop of size 4-8
    rseq_hairpin_matrix = [hairpin(rseq, i) for i in range(4,9)] 

    # matrix rows are the number of hairpin loop results (5 total for loop size 4,5,6,7,8)
    # matrix columns are the number of LSEQs / RSEQs -- each column is 1 sequence

    # flip the matrix on the diagonal so each sublist represents results for all loops of each sequence
    lseq_transposed_matrix = np.transpose(lseq_hairpin_matrix)
    rseq_transposed_matrix = np.transpose(rseq_hairpin_matrix)

    # transposed matrix rows/elements represent results for all loops of each sequence
    # transposed matrix columns represent the loop size - 4  --> to get loop size take column numbers + 4
    # need to compare each row (seq results) now, if any are True, that seq contains a hairpin with loop size(s): column size(s) + 4 

    # list stores loop sizes
    lseq_loop_sizes = [[i+4 for i,col in enumerate(sub_list) if col == True] for sub_list in lseq_transposed_matrix]
    rseq_loop_sizes = [[i+4 for i,col in enumerate(sub_list) if col == True] for sub_list in rseq_transposed_matrix]

    # finds max loop size of the present ones, if empty sublist that means there is no hairpin loop so the max size is 0
    lseq_max_loop = [max(loop_values) if len(loop_values) > 0 else 0 for loop_values in lseq_loop_sizes]
    rseq_max_loop = [max(loop_values) if len(loop_values) > 0 else 0 for loop_values in rseq_loop_sizes]


    # add lseq_max_loop to lseq df and rseq_max_loop to rseq df
    df_lseq = pd.concat([df_lseq, pd.DataFrame(lseq_max_loop, columns=['L_HAIRPIN'])], axis=1)
    df_rseq = pd.concat([df_rseq, pd.DataFrame(rseq_max_loop, columns=['R_HAIRPIN'])], axis=1)

    ##############################
    # drop columns that exist in df_lseq so joining the 2 dfs do not result in duplicate columns
    df_rseq = df_rseq.drop(columns=['CHROM', 'POS', 'REF', 'ALT'])

    # join df_lseq and df_rseq into a single df
    df = df_lseq.join(df_rseq)
    
    return df