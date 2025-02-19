#!/bin/bash
# activate a python 3.8 or 3.9 env with numpy and pandas installed
# example command: bash 10Xmapping_processing_to_matrix_v2.py matrix_table.tsv
# takes a table as input
# the input table has no header.
# the input file has the following columns
# sample_sn sample_wes path_of_parsed.txt path_of_maf path_of_Barcode_Annotation.txt

path=`pwd`
while read sample_sn sample_wes parsed maf barcodes; do
    mkdir -p ${path}/${sample_wes}
    echo ${sample_sn}
    echo ${sample_wes}
    echo ${parsed}
    echo ${maf}
    python3 /diskmnt/Projects/Users/austins2/tools/10Xmapping/10Xmapping_processing_to_matrix_v4.py ${barcodes} ${parsed} ${maf} 
done < $1
