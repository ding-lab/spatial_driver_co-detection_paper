#!/bin/bash
#example call:
#bash script.sh -t table.tsv -o /path/to/output/matrix/folder
#bash script.sh -t /diskmnt/Projects/snATAC_analysis/inferCNV-ATAC/CESC/outputs/matrices/CESC-infercnv-outs.txt -o .

while getopts "t:o:" opt; do
  case $opt in
    t)
      table=$OPTARG
      ;;
    o)
      outdir=$OPTARG
      ;;
    \?)
      echo "script usage: $(basename $0) [-t] " >&2
      exit 1
      ;;
  esac
done
while read sample observations tumor references; do
    bsub -G compute-dinglab -q general -J ${sample} -n 1 -oo ${outdir}/${sample}.map.log -R "select[mem>60GB] span[hosts=1] rusage[mem=60GB]" -M 60GB -a 'docker(austins2/condapython:3.9.6.base)' python3 /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/infercnv_postprocessing_v3.py ${observations} ${references} ${sample} ${outdir}
done < ${table}

