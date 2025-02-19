#!/bin/bash
# activate a python 3.8 or 3.9 env with numpy and pandas installed
# example command: bash 10Xmapping_processing_to_matrix_v2.py matrix_table.tsv
# takes a table as input
# the input table has no header.
# the input file has the following columns
# sample_sn sample_wes path_of_parsed.txt path_of_maf path_of_Barcode_Annotation.txt

path=`pwd`
while read sample maf var ref; do
    mkdir -p ${path}/${sample}
    python3 /diskmnt/Projects/Users/austins2/tools/10Xmapping/maf_from_matrix.py ${sample} ${maf} ${var} ${ref} "${path}/${sample}" 
done < $1
