#!/bin/bash
#example command: bash mapping_wrapper.sh mapping_table.tsv
#takes a table as input
#the input table has no header.
#the input file has the following columns
# snRNA_sample_ID WES_sample_ID bam_from_cellranger maf_from_somaticwrapper

path=`pwd`
#mkdir -p ${path}/${sample_wes}
while read sample_sn sample_wes bam maf; do
    mkdir -p ${path}/${sample_wes}  
    perl /diskmnt/Projects/Users/austins2/tools/10Xmapping/10Xmapping.pl --mapq 0 --bam ${bam} --maf ${maf} --out ${path}/${sample_wes}/${sample_sn}.mapped.txt
done < $1

