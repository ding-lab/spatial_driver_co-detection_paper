#!/bin/bash
while getopts "t:o:g:" opt; do
  case $opt in
    t )
      table="$OPTARG"
      ;;
    o )
      outdir="$OPTARG"
      ;;
    g )
      genes="$OPTARG"
      ;;
    \? )
      echo "script usage: $(basename $0) [-t] [-o] [-g]" >&2
      exit 1
      ;;
  esac
done
#mkdir -p ${outdir}/output/matrices
#source activate r4-sig

while read sample rds cnv_matrix_file; do
    echo ${sample}
    echo ${cnv_matrix_file}
    echo ${outdir}
    echo ${rds}
    echo ${genes}
    Rscript /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/plotting_cnvs_per_cell_type.R -i ${cnv_matrix_file} -s ${sample} -o ${outdir}  -r ${rds} -g ${genes}
done < ${table}
