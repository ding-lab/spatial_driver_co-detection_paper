#!/bin/bash
# example command: bash bamrc_wrapper_cp1.sh bam_table.tsv
# takes a table as input
# the input table has no header.
# the input file has the following columns
# snRNA_sample_ID WES_sample_ID path_to_tumor_wes_bam
# by: Austin Southard-smith
# based on: run_bamreadcount.sh by: Fernanda Martins Rodrigues (fernanda)
# last chagned: 20220615
while getopts "I:l:R:o:O:" opt; do
    case $opt in
        I)
            INPUT_BAM_TABLE=$OPTARG # input BAM file
            ;;
        l)
            SITE_LIST=$OPTARG # file containing a list of regions to report readcounts within
            ;;
        R)
            REF=$OPTARG # reference genome .fa file
            ;;
        O)
            OUT_DIR=$OPTARG # output directory
            ;;
        \?)
            echo "Invalid option -$OPTARG"
            echo $USAGE
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            echo $USAGE
            exit 1
            ;;
    esac
done

if [ ! -d ${OUT_DIR} ]; then
    mkdir ${OUT_DIR}
fi

mkdir ${OUT_DIR}/logs

#while read SAMPLE_SN SAMPLE_WES BAM; do
#    mkdir ${OUT_DIR}/${SAMPLE_WES}
#    mkdir ${OUT_DIR}/logs
#    LSF_DOCKER_ENTRYPOINT=/bin/bash bsub -G compute-dinglab -q dinglab -oo ${OUT_DIR}/logs/${SAMPLE_WES}.bamrc.log -a 'docker(mgibio/bam-readcount)' "bam-readcount -q 10 -b 15 -f ${REF} -l ${SITE_LIST} ${BAM} > ${OUT_DIR}/${SAMPLE_WES}/${SAMPLE_WES}.readcounts.tsv"
#done < ${INPUT_BAM_TABLE}

while read SAMPLE_SN SAMPLE_WES BAM; do
    mkdir ${OUT_DIR}/${SAMPLE_WES}
    STORAGE1="${OUT_DIR}/${SAMPLE_WES}"
    STORAGE2=`dirname ${BAM}`
    STORAGE3=`dirname ${REF}`
    STORAGE4=`dirname ${SITE_LIST}`
    LSF_DOCKER_VOLUMES="$HOME:$HOME ${STORAGE1}:${STORAGE1} ${STORAGE2}:${STORAGE2} ${STORAGE3}:${STORAGE3} ${STORAGE4}:${STORAGE4}" LSF_DOCKER_ENTRYPOINT=/bin/bash bsub -G compute-dinglab -q dinglab -oo ${OUT_DIR}/logs/${SAMPLE_WES}.bamrc.log -a 'docker(mgibio/bam-readcount)' "bam-readcount -q 10 -b 15 -f ${REF} -l ${SITE_LIST} ${BAM} > ${OUT_DIR}/${SAMPLE_WES}/${SAMPLE_WES}.readcounts.tsv"
done < ${INPUT_BAM_TABLE}
