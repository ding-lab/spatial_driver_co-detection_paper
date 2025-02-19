#!/bin/bash
#example command: bash parsing_wrapper.sh parsing_table.tsv
#takes a table as input
#the input table has no header.
#the input file has the following columns
# snRNA_sample_ID WES_sample_ID path_to_file_for_mutations

path=`pwd`
mkdir -p ${path}/${sample_wes}
while read sample_sn sample_wes mapping; do
    perl /diskmnt/Projects/Users/austins2/tools/10Xmapping/parse_scrna_bc.pl ${mapping} ${path}/${sample_wes}/${sample_sn}.parsed.txt
done < $1
