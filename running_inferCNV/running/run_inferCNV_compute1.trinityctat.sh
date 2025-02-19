#!/bin/bash
#example call:
#bash run_inferCNV_compute1.sh -T reference_annotation.tsv -D run_dir -G job_group 
#bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/run_inferCNV_compute1.trinityctat.sh -T /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-06-06/inputs/annotations_file/reference_cells.txt -D /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-06-06 -G /a.n.southard-smith/infercnv


while getopts "T:D:G:" opt; do
  case $opt in
    T)
      table=$OPTARG
      ;;
    D)
      rundir=$OPTARG
      ;;
    G)
      jgroup=$OPTARG
      ;;
    \?)
      echo "script usage: $(basename $0) [-T] [-D] " >&2
      exit 1
      ;;
  esac
done

while read sample reference tumor; do
    dir_run=${rundir}
    ## specify the docker image, this is an environment with all the dependent R packages installed
    docker_image_name=trinityctat/infercnv
    ## specify the inferCNV running mode
    analysis_mode=subclusters
    ## specify the path to the raw count matrix
    ### specify the aliquot id to enable getting the path to the raw count matrix
    ### CHANGE THIS for each aliquot
    snRNA_aliquot_id=${sample}
    ### specify the directory with all the raw count matrix
    dir_raw_count_file=${rundir}/inputs/raw_counts_matrix/
    ### specify the path to the raw count matrix for this aliquot
    path_raw_count_file=${dir_raw_count_file}${snRNA_aliquot_id}.raw_counts.tsv
    ## specify the path to the annotation file
    ### specify the run id for this run because the annotation file might change in the future
    now=`date +'%Y%m%d'`
    run_id=${now}
    ### specify the directory with all the annotation files
    dir_annotation_file=${rundir}/inputs/annotations_file/${snRNA_aliquot_id}/
    ### specify the path with the annotation file for this aliquot
    path_annotation_file=${dir_annotation_file}${snRNA_aliquot_id}.Barcode_Annotation.txt
    ## specify the path to the gene order file, this file contains the position information of each gene
    path_gene_order_file=/storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/genes.Cellranger-2020-A.refdata-gex-GRCh38-2020-A.gene_filterd.remove_excess_columns.trimmed.no_duplicates.txt
    ## specify the cutoff for the average gene “expression” level to filter genes
    cutoff=0.1 #this was originally set a 0.04 from Nadya
    ## specify the directory to the outputs for this aliquot
    ### specify the directory to the outputs for all the runs
    dir_infercnv_outputs=${rundir}/outputs/
    ### specify the directory to the outputs for this run
    dir_infercnv_outputs_by_run=${dir_infercnv_outputs}${run_id}/
    mkdir -p ${dir_infercnv_outputs_by_run}
    ### specify the directory to the outputs for this run and for this aliquot
    dir_output=${dir_infercnv_outputs_by_run}${snRNA_aliquot_id}/
    mkdir -p ${dir_output}
    ## specify the label for the reference cells in the annotation file
    ### for reference cells with different cell types, please combine them here by “,”, for example “Endothelial cells,Macrophages,CD4+ T-cells"
    ### if grouping all of the reference cells together under a single "normal" lable then proved that here as just a single label. example: "Normal"
    ref_group_names=${reference}
    ## specify the path to the log file, each run has a unique timestamp
    path_log_file=${dir_output}${snRNA_aliquot_id}.$(date +%Y%m%d%H%M%S).log
    ## below is the main running code
    ### the following line specify 
    echo ${ref_group_names}
    export LSF_DOCKER_PRESERVE_ENVIRONMENT=false 
    bsub -G compute-dinglab -q dinglab -g ${jgroup} -J ${tumor} -N -n 12 -R 'select[mem>240GB] rusage[mem=240GB]' -M 240GB -oo ${path_log_file} -a 'docker(trinityctat/infercnv:1.11.2)' "ulimit -s unlimited && Rscript /storage1/fs1/dinglab/Active/Projects/austins2/pancan-ATAC/matching-snRNAseq/PDAC/inferCNV/inferCNV.trinityctat.R \
            --analysis_mode=${analysis_mode} \
            --raw_counts_matrix=${path_raw_count_file} \
            --annotations_file=${path_annotation_file} \
            --gene_order_file=${path_gene_order_file} \
            --cutoff=${cutoff} \
            --out_dir=${dir_output} \
            --cluster_by_groups \
            --denoise \
            --HMM \
            --num_threads=12 \
            --ref_group_names=${ref_group_names}"
done < ${table}

