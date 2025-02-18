#!/usr/bin/env python3
#By: Austin Southard-Smith
# example command:
# conda activate py3.9
# python3 bamrc_to_vaf_table_v3.py \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/variant_look_up_table.tsv \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/input_bamrc_counts_files.tsv \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/ \

import numpy as np
import pandas as pd
import sys
import os

path_to_variant_translation_table = sys.argv[1]
path_to_bamrc_input_table = sys.argv[2] # this is the list of samples that will be included in the output table.
out_dir = sys.argv[3]

# read in a table of variants can just use the translation layer
# index	Chromosome	Start	Stop	Ref_allele	Alt_allele	HGVSc	hugo_symbol	HGVSp_short	hugo_HGVSp	bamrc_parsed_alt_allele	bamrc_parsed_ref_allele	alt_probe_feature	ref_probe_feature	matrix_key	Chromosome__Start_key
# can just read this in a a pd df
# input should have a colnames and index column of numbers
# only include variants that will be present in the output table.
variant_table = pd.read_csv(path_to_variant_translation_table, sep = "\t", low_memory=False, header=0, index_col=0)
print(variant_table.columns)
print(variant_table.index)
hugo_HGVSp_list = variant_table['hugo_HGVSp'].tolist()
chrom_start_key_list = variant_table["Chromosome__Start_key"].tolist()
# read in a table of samples files
# the output should be a dictionary {"sample": "output_folder"}
# input should have a colnames and index column of numbers
# index	Truncated_Xenium_ID	Xenium_ID	bulk_WES_ID_tumor	bulk_WES_ID_normal	bamrc_tumor_compute1_path	bamrc_normal_compute1_path	bamrc_tumor_katmai_path	bamrc_normal_katmai_path
# only include samples that will be present in the output table.
sample_table = pd.read_csv(path_to_bamrc_input_table, sep = "\t", low_memory=False, header=0, index_col=0)

sample_list = sample_table["Xenium_ID"].tolist()

# make the vaf table
vaf_df = pd.DataFrame(np.nan, index=sample_list, columns=hugo_HGVSp_list, dtype = np.float64)

print(vaf_df)
complete_list_of_bamrc_lines = []
max_line_length = 0
# function to parse the parsed bamrc files
def bamrc_parse(df = vaf_df, path_to_bamrc = "/path/to/bamrc", var_table = pd.DataFrame(), samp = "sample_1",list_of_bamrc_lines = [], line_max = 0):
    bamrc_file = open(file = path_to_bamrc, mode = 'r') #open the file
    #there is no header
    is_pdac = False
    if samp in ["HT227P1-S1H1L1U1","HT242P1-S1H4L4U1","HT284P1-S1H1A1U1","SP001P1-Fp1U1","HT270P1-S1H1A1US2_1","HT060P1-S1R1Fp1U1","HT061P1-S1P1A1L1U1","HT061P1-S1P1A1L4U1","HT125P1-S1H4A1L1U1","HT125P1-S1H8A1U1","HT185P1-S1H2L1U1"]:
        is_pdac = True
    for line in bamrc_file:
        line_list = line.strip().split()
        line_list.insert(0,samp)
        if len(line_list) > line_max:
            line_max = len(line_list)
        list_of_bamrc_lines.append(line_list)
        chrom_start_key = line_list[1] + "__" + str(line_list[2])
        if chrom_start_key in chrom_start_key_list:
            print(chrom_start_key)
            coverage = int(line_list[4])
            allele_dict = dict()
            for i in range(5, len(line_list)-1, 2):
                allele = line_list[i]
                allele_count = int(line_list[i+1])
                # print(allele)
                # print(allele_count)
                allele_dict[allele] = allele_count
            bamrc_parsed_alt_allele_list = list(var_table.loc[(var_table["Chromosome__Start_key"] == chrom_start_key),"bamrc_parsed_alt_allele"])
            print(bamrc_parsed_alt_allele_list)
            for bamrc_parsed_alt_allele in bamrc_parsed_alt_allele_list:
                if bamrc_parsed_alt_allele in list(allele_dict.keys()):
                    alt_allele_count = allele_dict[bamrc_parsed_alt_allele]
                    hugo_HGVSp_key = var_table.loc[(var_table["Chromosome__Start_key"] == chrom_start_key) & (var_table["bamrc_parsed_alt_allele"] == bamrc_parsed_alt_allele), "hugo_HGVSp"].item()
                    # print(hugo_HGVSp_key)
                    # print(sample)
                    # print(alt_allele_count)
                    # print(coverage)
                    if is_pdac: 
                        if alt_allele_count > 1:
                            df.at[samp, hugo_HGVSp_key] = alt_allele_count/coverage
                        else:
                            df.at[samp, hugo_HGVSp_key] = 0/coverage
                    else:
                        if alt_allele_count > 3:
                            df.at[samp, hugo_HGVSp_key] = alt_allele_count/coverage
                        else:
                            df.at[samp, hugo_HGVSp_key] = 0/coverage
                else:
                    alt_allele_count = 0 #this handles insertions and deletions if they are not present in the bamrc output.
                    hugo_HGVSp_key = var_table.loc[(var_table["Chromosome__Start_key"] == chrom_start_key), "hugo_HGVSp"].item()
                    # print(hugo_HGVSp_key)
                    # print(sample)
                    # print(alt_allele_count)
                    # print(coverage)
                    df.at[samp, hugo_HGVSp_key] = alt_allele_count/coverage
            # need to account for insertions and deletions
            # positive strand deletion will be listed at the start position as "-[ref_allele]"
            # negative strand deletion will be listed at the start position as "-[ref_allele]"
            # positive strand insertion will be listed at the start position as "+[Alt_allele]"
            # we didn't test any negative strand insertions
    return df.copy(), list_of_bamrc_lines, line_max

#read in list of bamrc files
#loop over each file and add contents of it to table.
for sample in sample_list:
    print(sample)
    bamrc_output_path = str(sample_table.loc[(sample_table["Xenium_ID"] == sample),"bamrc_tumor_katmai_path"].item())
    print(bamrc_output_path)
    if not sample_table.loc[(sample_table["Xenium_ID"] == sample),"bamrc_tumor_katmai_path"].isnull().item():
        vaf_df, complete_list_of_bamrc_lines, max_line_length = bamrc_parse(df = vaf_df, path_to_bamrc = bamrc_output_path, var_table = variant_table, samp = sample, line_max = max_line_length)

vaf_df.insert(loc=0, column="index", value=sample_list)
print(vaf_df)
vaf_df.to_csv(out_dir + "/" +"bamrc_vaf_summary_table.tsv", sep = "\t" ,index=False, header=True)

bamrc_out_file = open(out_dir + "/" +"complete_parsed_bamrc_result_file.tsv", mode = 'w')
bamrc_out_file.write("sample"+'\t'+'chromosome'+'\t'+'start'+'\t'+'reference_allele'+'\t'+'coverage'+'\t'+'A'+'\t'+'counts_of_A'+'\t'+'C'+'\t'+'counts_of_C'+'\t'+'G'+'\t'+'counts_of_G'+'\t'+'T'+'\t'+'counts_of_T'+'\t'+'indel_if_present'+'\t'+'counts_of_indel'+'\n')
"1"+'\t'+'2'+'\t'+'3'+'\t'+'4'+'\t'+'5'+'\t'+'6'+'\t'+'7'+'8'+'\t'+'9'+'10'+'\t'+'11'+'12'+'\t'+'13'+'14'+'\t'+'15'
["sample",'chromosome','start','reference_allele','coverage','A','counts_of_A','C','counts_of_C','G','counts_of_G','T','counts_of_T','indel_if_present','counts_of_indel']
bamrc_out_file.writelines("")
for line in complete_list_of_bamrc_lines:
    new_line = [] 
    for i in range(max_line_length):
        if i > (len(line)-1):
            value_to_append = "NA"
        else:
            value_to_append = line[i]
        if (i == max_line_length-1):
            new_line.append(value_to_append+'\n')
        else:
            new_line.append(value_to_append+'\t')
    bamrc_out_file.writelines(new_line)
bamrc_out_file.close()
