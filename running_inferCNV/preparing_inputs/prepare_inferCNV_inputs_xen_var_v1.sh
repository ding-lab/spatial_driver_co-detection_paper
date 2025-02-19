#!/bin/bash
#example call:
#bash script.sh -t sample_table.tsv -o ouput_dir -m seurat_meta.data_table.tsv

while getopts "t:o:m:" opt; do
  case $opt in
    t)
      table=$OPTARG
      ;;
    o)
      outdir=$OPTARG
      ;;
    m)
      metadata=$OPTARG
      ;;
    \?)
      echo "script usage: $(basename $0) [-t] [-o] " >&2
      exit 1
      ;;
  esac
done
mkdir -p ${outdir}/inputs
mkdir -p ${outdir}/outputs
while read sample rds_file tumor; do
    Rscript /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/prepare_inferCNV_inputs_xen_var_v1.R -i ${rds_file} -s ${sample} -o "${outdir}/inputs" -t ${tumor} #-m ${metadata} 
done < ${table}

