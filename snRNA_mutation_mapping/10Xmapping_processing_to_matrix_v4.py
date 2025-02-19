#!/usr/bin/env python3
'''
By: Austin Southard-Smith
usage: python3 10Xmapping_processing_to_matrix.py barcode_annotations.txt  mutations.parsed.txt sample_maf_file.tsv
Example:

python3 10Xmapping_processing_to_matrix.py /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/BRCA/inferCNV/inputs/annotations_file/HT206B1-XBN1-1/HT206B1-XBN1-1.Barcode_Annotation.txt mutations_10Xmapped_HT206B1-XBN1-1.parsed.txt /diskmnt/Projects/ATAC_pancan_analysis/bulk_somatic_mutation/BRCA/HT206B1-TWCE-HT206B1-S1H4A2Y2D1_1_T.dnp.annotated.maf /path/to/output/dir


argv[1] <-- barcode_annotations.txt = tsv file with cell barcodes in column one and anything else in any other columns (we only care about column 1 here) (this file is also used as input for inferCNV)
argv[2] <-- mutations.parsed.txt = output of the second step of running the 10Xmapping pipeline
looks like this:
chrChromosome   Start_Position  Reference_Allele        Tumor_Seq_Allele2       Ref-Support
chrChromosome   Start_Position  Reference_Allele        Tumor_Seq_Allele2       Var-Support
chr1    1922998 G       A       Ref-Support
chr1    1922998 G       A       Var-Support <- variant allele
chr1    6461492 C       T       Ref-Support <- reference allele
A00325:357:HTN5KDSXY:4:1102:13747:1908  AGTACCAAGAGTAACT-1 <- barcodes with reference allele
A00325:357:HTN5KDSXY:4:1102:14073:1971  AGTACCAAGAGTAACT-1

argv[3] <-- samples specific maf file
'''

import os
import sys
import numpy as np
import pandas as pd

#constructing the output filename:
print(0,sys.argv[0])
print(1,sys.argv[1])
print(2,sys.argv[2])
print(3,sys.argv[3])
basename = os.path.splitext(os.path.basename(sys.argv[2]))[0]
output_variant_file_name=sys.argv[3]+"/"+basename+".mutation_matrix.tsv"
output_variant_gene_file_name=sys.argv[3]+"/"+basename+".variant_gene_matrix.tsv"
output_reference_gene_file_name=sys.argv[3]+"/"+basename+".reference_gene_matrix.tsv"

sample_maf_file_path = "path/to/sample/sample.Barcode_Annotation.txt" #this is the same barcode annotation file that is used for inferCNV
sample_barcode_path = sys.argv[1]
sample_barcodes_df = pd.read_csv(sample_barcode_path, sep="\t", header = None, index_col = None) #reading in the same barcode annotation file used for inferCNV (see output of /diskmnt/Projects/Users/austins2/tools/inferCNV/prepare_inferCNV_inputs.R) -- column1 is cell barcodes and column 2 is cell-type annotations
sample_barcodes_df = sample_barcodes_df.transpose() #now the first row is cell barcodes and 2nd is cell types
sample_barcodes_df.columns = sample_barcodes_df.iloc[0,:] #now the cell barcodes is the header (Cell barcodes are still also in row 1
sample_barcodes_df = sample_barcodes_df.drop(sample_barcodes_df.index[1]) #removing the cell types in the 2nd row (1 based)
#sample_barcodes_df = sample_barcodes_df.drop(sample_barcodes_df.index[0]) #this line will drop the first row which is sample barcodes. <- don't do this yet since then the dataframe would be empty columns.
whitelist_barcodes = list(sample_barcodes_df.columns)
#read your maf file into a dictionary that will include the variant code and the position/description of the variant {Chromosome_Start_Position_Reference-Allele_Tumor-Seq-Allele2: Gene_HGVSp}
sample_maf_file_path = "path/to/sample/maf_file" #/diskmnt/Projects/ATAC_pancan_analysis/bulk_somatic_mutation/BRCA/HT206B1-TWCE-HT206B1-S1H4A2Y2D1_1_T.dnp.annotated.maf
sample_maf_file_path = sys.argv[3]
maf_file = open(sample_maf_file_path)
header = maf_file.readline()
maf_dict = {}
gene_dict = {}
gene_set = set()
maf_dict_detailed = {}
for line in maf_file:
    line_list = line.strip().split("\t")
    variant_key = "_".join([line_list[4], line_list[5], line_list[10], line_list[12]]) #Chromosome_Start_Position_Reference-Allele_Tumor-Seq-Allele2 is the key from the maf file
    variant_code = "_____".join([line_list[0], line_list[36]]) #Gene_HGVSp is the code from the maf file where Gene is the column referred to as Hugo_Symbol. (Note: due to base wobble it is possible for multiple variants to result in the same variant code if only using Gene___HGVSp. In order to truly see each genetic variant split out into it's own column then you will need to change this line so that it results only in unique variant codes. An example of one that would do this is variant_code = "_____".join([line_list[0], line_list[36], line_list[4], line_list[5], line_list[10], line_list[12]]) . This results in a line formatted like "KRAS_____p.G12D_____chr12_____25245350_____25245350_____C_____T" . It is necessary to format the key like this due to the multiplicity of symbols used in both Gene symbols and HGVSp. This ensures future analysis can recover each unique portion of the variant code is a simple manner
    variant_code_detailed = "_____".join([line_list[0], line_list[36], line_list[4], line_list[5], line_list[10], line_list[12]])
    maf_dict[variant_key] = variant_code  #{Chromosome_Start_Position_Reference-Allele_Tumor-Seq-Allele2: Gene_HGVSp} <- the separator here in the code is 5 underscores in a row

    maf_dict_detailed[variant_key] = variant_code_detailed
    gene_dict[variant_key] = line_list[0]
    gene_set.add(line_list[0])
maf_file.close()


#now expand your initial data sample_barcodes dataframe so that it has 1 row for every variant.
num_variants = len(list(maf_dict.values()))
num_barcodes = len(sample_barcodes_df.iloc[0])
num_genes = len(list(gene_set))
num_variants_detailed = len(list(maf_dict_detailed.values()))
print(num_genes)

def add_zeros_array_to_pd_dataframe_as_new_rows(df,length, width):
    array_to_add_to_barcodes = np.zeros((length,width)) #setup the 0s array
    new_df = pd.DataFrame(array_to_add_to_barcodes, index = df.columns) #convert to pd dataframe
    new_df = new_df.astype(int)
    return new_df

barcode_variant_df = add_zeros_array_to_pd_dataframe_as_new_rows(df = sample_barcodes_df, length = num_barcodes, width = num_variants)

barcode_gene_reference_allele_df = add_zeros_array_to_pd_dataframe_as_new_rows(df = sample_barcodes_df, length = num_barcodes, width = num_genes)
barcode_gene_variant_allele_df = add_zeros_array_to_pd_dataframe_as_new_rows(df = sample_barcodes_df, length = num_barcodes, width = num_genes)
barcode_gene_mutation_df = add_zeros_array_to_pd_dataframe_as_new_rows(df = sample_barcodes_df, length = num_barcodes, width = num_genes)


barcode_variant_allele_df_detailed = add_zeros_array_to_pd_dataframe_as_new_rows(df = sample_barcodes_df, length = num_barcodes, width = num_variants_detailed)
barcode_reference_allele_df_detailed = add_zeros_array_to_pd_dataframe_as_new_rows(df = sample_barcodes_df, length = num_barcodes, width = num_variants_detailed)


#assign variant codes as row names in barcode_variant dataframe
var_codes = sorted(list(maf_dict.values()))
barcode_variant_df.columns = var_codes
barcode_gene_reference_allele_df.columns = sorted(list(gene_set))
barcode_gene_variant_allele_df.columns = sorted(list(gene_set))
barcode_gene_mutation_df.columns = sorted(list(gene_set))
var_codes_detailed = sorted(list(maf_dict_detailed.values()))
barcode_variant_allele_df_detailed.columns = var_codes_detailed
barcode_reference_allele_df_detailed.columns = var_codes_detailed


#read your variants into a dictionary
sample_path = "/path/to/sample/sample.mutations.parsed.txt"
sample_path = sys.argv[2]
sample_id = sample_path.split('/')
variant_dictionary={} #this will be a nested dictionary: {Variant_key:{"Reference":[barcode1,barcode2,barcode3],"Variant":[barcode1,barcode2,barcode3]}}
gene_barcode_dictionary = {} #this will be a nested dictionary: {Gene_key:{"Reference":[barcode1,barcode2,barcode3],"Variant":[barcode1,barcode2,barcode3]}}
current_variant = ""
current_key = "" #this will be the key (either 'Reference' or Variant)
current_gene = "" 
sample_file = open(sample_path, 'r')
for line in sample_file:
    if "Chromosome" in line:
        continue
    elif "Support" in line:
        chrom,start,ref_allele,var_allele,allele_support = line.strip().split()
        variant = chrom + "_" + start + "_" + ref_allele + "_" + var_allele
        gene = gene_dict[variant]
        if allele_support == "Ref-Support":
            current_key = "Reference" #this will 
        else:
            current_key = "Variant"
        if current_variant != variant:
            current_variant = variant
            variant_dictionary[variant] = {"Reference":[], "Variant":[]}
        if current_gene != gene:
            current_gene = gene
            gene_barcode_dictionary[gene] = {"Reference":[], "Variant":[]}
    else:
        flow_cell_info,cell_barcode = line.strip().split()
        if cell_barcode in whitelist_barcodes:
            gene_barcode_dictionary[gene][current_key].append(cell_barcode)
            variant_dictionary[current_variant][current_key].append(cell_barcode)

sample_file.close()

#print("gene_barcode_dictionary")
#print(gene_barcode_dictionary)

#print("variants")
#print(sorted(list(variant_dictionary.keys())))

#loop over the varaint_dictionary and assign the positions at different barcodes based on if the alternate allele is present or not.
#0 = unknown and variant site not sequenced in data
#1 = reference allele present in sequence data
#2 = variant allele present in sequence data


variants = sorted(list(variant_dictionary.keys()))
for var in variants:
    gene = gene_dict[var]
    print(gene)
    for barcode in variant_dictionary[var]["Reference"]:
        barcode_variant_df.at[barcode, maf_dict[var]] = 1
        barcode_variant_df = barcode_variant_df.copy()
    for barcode in variant_dictionary[var]["Variant"]:
        barcode_variant_df.at[barcode, maf_dict[var]] = 2
        barcode_variant_df = barcode_variant_df.copy()
    for barcode in gene_barcode_dictionary[gene]["Reference"]:
        barcode_gene_mutation_df.at[barcode,gene] = 1
        barcode_gene_mutation_df.copy()
    for barcode in gene_barcode_dictionary[gene]["Variant"]:
        barcode_gene_mutation_df.at[barcode,gene] = 2
        barcode_gene_mutation_df.copy()
for gene in sorted(list(gene_barcode_dictionary.keys())):
    variant_barcodes = sorted(list(gene_barcode_dictionary[gene]["Variant"]))
    # print(variant_barcodes)
    reference_barcodes = sorted(list(gene_barcode_dictionary[gene]["Reference"]))
    # print(reference_barcodes)
    for barcode in variant_barcodes:
        # print(gene, barcode)
        barcode_gene_variant_allele_df.at[barcode,gene] += 1
        barcode_gene_variant_allele_df.copy()
        # print(barcode_gene_variant_allele_df.at[barcode,gene])
    for barcode in reference_barcodes:
        # print(gene, barcode)
        barcode_gene_reference_allele_df.at[barcode,gene] += 1
        barcode_gene_reference_allele_df.copy()
for var in variants:
    variant_barcodes = sorted(list(variant_dictionary[var]["Variant"]))
    # print(variant_barcodes)
    reference_barcodes = sorted(list(variant_dictionary[var]["Reference"]))
    for barcode in variant_barcodes:
        barcode_variant_allele_df_detailed.at[barcode, maf_dict_detailed[var]] += 1
        barcode_variant_allele_df_detailed.copy()
    for barcode in reference_barcodes:
        barcode_reference_allele_df_detailed.at[barcode, maf_dict_detailed[var]] += 1
        barcode_reference_allele_df_detailed.copy()


#print(barcode_variant_df)
#print(barcode_gene_reference_allele_df)
#print(barcode_gene_variant_allele_df)

#write the final dataframe to a tsv file that is compatible with R and can be added as metadata in Seurat
current_working_dir = os.getcwd()

# In this output matrix a 0 is no allele (reference or variant detected), 1 is reference allele detected, 2 is variant allele detected. This matrix reports one column per gene per unique HGVSp present in the maf file.
# If there are multiple variants present that can result in the same Gene___HGVSp key (such as may result from base wobble) then this column will combine the results of all of them to this column. To change this behavior simply go up to where the maf dictionary is constructed (lines 54-61) and change how the variant_code is made on line 57 by choosing more or different values from the maf file to define a variant. Alternatively, un-comment line 58 and comment out line 57. 
output_variant_file_name=current_working_dir+"/"+basename+".mutant_variant_matrix.tsv"
# In this output matrix a 0 is no allele (reference or variant detected), 1 is reference allele detected, 2 is variant allele detected. This matrix only reports one column per gene (Hugo_Symbol) from the maf file.
output_mutant_gene_file_name=current_working_dir+"/"+basename+".mutation_matrix.tsv"
# In this output matrix the number of variants alleles for a gene are reported. If multiple variant sites or variants at the same site are reported for each gene they are not split out here.
output_variant_gene_file_name=current_working_dir+"/"+basename+".variant_gene_matrix.tsv"
# In this output matrix the number of reference alleles for a gene are reported. If multiple variant sites or variants at the same site are reported for each gene they are not split out here.
output_reference_gene_file_name=current_working_dir+"/"+basename+".reference_gene_matrix.tsv"

output_variant_allele_count_file_name=current_working_dir+"/"+basename+".variant_allele_count_matrix.tsv"
output_reference_allele_count_file_name=current_working_dir+"/"+basename+".reference_allele_count_matrix.tsv"


#output = barcode_variant_df.transpose(copy=True)
#output.columns = output.iloc[0]
output = barcode_variant_df.copy()
#output_variant_gene = barcode_gene_variant_allele_df.transpose(copy=True)
#output_variant_gene.columns = output_variant_gene.iloc[0]
output_reference_gene = barcode_gene_reference_allele_df.copy()
#output_reference_gene = barcode_gene_reference_allele_df.transpose(copy=True)
#output_reference_gene = output_reference_gene.iloc[0]
output_variant_gene = barcode_gene_variant_allele_df.copy()

output_mutant_gene = barcode_gene_mutation_df.copy()

output_variant_allele_count = barcode_variant_allele_df_detailed.copy()
output_reference_allele_count = barcode_reference_allele_df_detailed.copy()

#output = output.drop(output.index[0])
# print(output)
# If there are multiple variants present that can result in the same Gene___HGVSp key (such as may result from base wobble) then this column will combine the results of all of them to this column. To change this behavior simply go up to where the maf dictionary is constructed (lines 54-61) and change how the variant_code is made on line 57 by choosing more or different values from the maf file to define a variant. Alternatively, un-comment line 58 and comment out line 57. 
output.to_csv(output_variant_file_name,sep = "\t" ,index=True, header=True)

# print(output_mutant_gene)
# In this output matrix a 0 is no allele (reference or variant detected), 1 is reference allele detected, 2 is variant allele detected. This matrix only reports one column per gene (Hugo_Symbol) from the maf file.
output_mutant_gene.to_csv(output_mutant_gene_file_name, sep = '\t', index = True, header = True)

# print(output_reference_gene)
# In this output matrix the number of variants alleles for a gene are reported. If multiple variant sites or variants at the same site are reported for each gene they are not split out here.
output_reference_gene.to_csv(output_reference_gene_file_name, sep='\t', index = True, header = True)

# print(output_variant_gene)
# In this output matrix the number of reference alleles for a gene are reported. If multiple variant sites or variants at the same site are reported for each gene they are not split out here.
output_variant_gene.to_csv(output_variant_gene_file_name, sep='\t',index = True, header=True)

 
output_variant_allele_count.to_csv(output_variant_allele_count_file_name, sep="\t", index = True, header = True)
output_reference_allele_count.to_csv(output_reference_allele_count_file_name, sep="\t", index = True, header = True)






