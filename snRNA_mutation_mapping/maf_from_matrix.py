import os
import sys
import pandas as pd

print(0, sys.argv[1]) # sample
print(1, sys.argv[2]) # maf
print(2, sys.argv[3]) # variant count table

print(3, sys.argv[4]) # reference count table
print(4, sys.argv[5]) # output_dir


basename = sys.argv[1]
out_dir = sys.argv[5]

path_ref_df = sys.argv[4]
ref_df = pd.read_csv(path_ref_df, sep='\t', low_memory=False, header=0, index_col="0")
path_var_df = sys.argv[3]
var_df = pd.read_csv(path_var_df, sep = '\t', low_memory=False, header=0,index_col="0")
path_maf = sys.argv[2]
maf_df = pd.read_csv(path_maf, sep = '\t', low_memory=False, header=0, index_col=None)


maf_indices = maf_df.index
#"KRAS_____p.G12D_____chr12_____25245350_____25245350_____C_____T" 


for idx in maf_indices:
    chromosome = maf_df.iloc[idx,4]
    start_position = maf_df.iloc[idx,5]
    ref_allele = maf_df.iloc[idx,10]
    alt_allele = maf_df.iloc[idx,12]
    gene = maf_df.iloc[idx,0]
    HGVSp_code = maf_df.iloc[idx,36]
    variant_code = "_____".join([str(gene), str(HGVSp_code), str(chromosome), str(start_position), str(ref_allele), str(alt_allele)])
    alt_count = var_df[variant_code].sum() # t_alt_count
    ref_count = ref_df[variant_code].sum() # t_ref_count
    maf_df.iloc[idx,39] = alt_count + ref_count # t_depth
    maf_df.iloc[idx,40] = ref_count
    maf_df.iloc[idx,41] = alt_count
    maf_df.iloc[idx,15] = basename
    maf_df.iloc[idx,16] = basename
    maf_df.copy()


output_allele_summary_file = out_dir+"/"+basename+".all_alleles_summary_output.maf"
maf_df.to_csv(output_allele_summary_file, sep='\t',index = False, header=True)


if int(maf_df["t_alt_count"].sum()) > 0:
    mutated_only_maf_df = maf_df.loc[(maf_df["t_alt_count"] != 0)]
    mutated_only_maf_df = mutated_only_maf_df.copy()
    output_allele_mutated_file = out_dir+"/"+basename+".mutated_only_output.maf"
    mutated_only_maf_df.to_csv(output_allele_mutated_file, sep='\t',index = False, header=True)



