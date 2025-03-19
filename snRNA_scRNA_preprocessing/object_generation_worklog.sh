# worklog for alignment and doublet detection via scrublet in previously unpublished datasets
# the different R and bash scripts are 
# running cellranger v6.0.2
export PATH=/diskmnt/Projects/Users/austins2/software/cellranger-6.0.2/bin:$PATH
# running cellranger-arc v2.0.0
export PATH=/diskmnt/Software/cellranger-arc-2.0.0/bin:$PATH
# HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1
cd /diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/20231030/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1/
cellranger-arc count --id=HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1 --disable-ui --localmem=180 --localcores=18 \
--libraries /diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/preprocessing/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1/TWCE-HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1_library.csv \
--reference /diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
# HT425B1-S1H1A3Y1N1Z1_1Bmn1_1
cd /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/
cellranger count --id HT425B1-S1H1A3Y1N1Z1_1Bmn1_1 --chemistry ARC-v1 --localcores 10 --localmem 50 --include-introns \
--transcriptome /diskmnt/Datasets/Reference/Cellranger-2020-A/refdata-gex-GRCh38-2020-A \
--fastqs ./preprocessing/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1/
# HT179C1-T1A3Y1N1
cd /diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/HT179C1-T1A3Y1N1/
cellranger-arc count --id HT179C1-T1A3Y1N1 --disable-ui --localcores 30 --localmem 200 \
--reference /diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--libraries /diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/preprocessing/HT179C1-T1A3Y1N1/library.csv
# C3L-03372_CPT0275960013_2022-01-31
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/cellranger-arc/
cellranger-arc count --id C3L-03372_CPT0275960013_2022-01-31 --disable-ui --localcores 32 --localmem 300 \
--reference /diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--libraries /diskmnt/Projects/GBM_sc_analysis/multiome_autoprocess/inputs/C3L-03372_CPT0275960013_2022-01-31/library.csv 
# C3N-00663_CPT0087730015_2021-03-08
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/cellranger-arc/
cellranger-arc count --id C3N-00663_CPT0087730015_2021-03-08 --disable-ui --localcores 32 --localmem 300 \
--reference /diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--libraries /diskmnt/Projects/GBM_sc_analysis/multiome_autoprocess/inputs/C3N-00663_CPT0087730015_2021-03-08/library.csv

# running scrublet worklog
# the scripts for running the available from the github repo https://github.com/Aust1nS2/automated_scrublet
# the environment used was the py3.9 environment from the github repo https://github.com/Aust1nS2/automated_scrublet
# HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030
mdkir ATAC combined RNA
# RNA portion can be run at the same time as the ATAC
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-RNA-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/20231030 \
/diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030/RNA
# ATAC portion can be run at the same time as the RNA
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-ATAC-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/20231030 \
/diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030/ATAC
# After ATAC and RNA are done then combine the results
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-auto-combining.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/20231030 \
/diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030/RNA \
/diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030/ATAC \
/diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030/combined
# HT425B1-S1H1A3Y1N1Z1_1Bmn1_1
cd /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/scrubblet/
mdkir RNA
cd /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/scrubblet/RNA
# Cellranger so RNA only
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-RNA-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_BRCA_primary/Cellranger \
/diskmnt/Projects/HTAN_BRCA_primary/Cellranger/scrubblet/RNA
# HT179C1-T1A3Y1N1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130/
mkdir ATAC combo RNA
# RNA portion can be run at the same time as the ATAC
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-RNA-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/ \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130//RNA
# ATAC portion can be run at the same time as the RNA
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-ATAC-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/ \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130/ATAC
# After ATAC and RNA are done then combine the results
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-auto-combining.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/ \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130/RNA \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130/ATAC \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130/combined
# C3L-03372_CPT0275960013_2022-01-31 and C3N-00663_CPT0087730015_2021-03-08
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/
mkdir ATAC combo RNA
# RNA portion can be run at the same time as the ATAC
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-RNA-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/cellranger-arc/ \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/RNA
# ATAC portion can be run at the same time as the RNA
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-ATAC-auto.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/cellranger-arc/ \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/ATAC
# After ATAC and RNA are done then combine the results
bash /diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts/scrublet-auto-combining.sh \
/diskmnt/Projects/Users/austins2/tools/automated_scrublet/scripts \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/cellranger-arc/ \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/RNA \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/ATAC \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/combined




# Individual Seurat objects were initially generated using the following scripts in seurat4 for separate projects. These objects were then uniformly re-processed in the Seurat 5.
# below is the command for the seurat 4 individual sample level object generation
# scripts and supporting files used here are present in `spatial_driver_co-detection_paper/snRNA_scRNA_preprocessing/scripts`
# HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1
conda activate r4.2 # this is the environment contained in `spatial_driver_co-detection_paper/snRNA_scRNA_preprocessing/snRNA_processing_worklog.sh`
cd /diskmnt/Projects/Users/austins2/pdac/primary_tumor/sn/individual-snRNA/SCTv2_with_doublets
bash /diskmnt/Projects/Users/austins2/tools/seurat4-scrublet-no-filter-sctv2-cellranger-arc2.sh /diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/20231030/ /diskmnt/Projects/HTAN_analysis_2/PDAC/primary/sn/scrublet-arc/20231030/combined
# HT425B1-S1H1A3Y1N1Z1_1Bmn1_1
conda activate r4.2
/diskmnt/Projects/HTAN_BRCA_primary/Cellranger/logs/seurat4_no_removal_scrublet_cellranger6_multiomeRNA.out
# /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/seurat/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1_processed.rds
# Rscript /diskmnt/Projects/Users/fernanda/Scripts/Cellranger/PowerRanger/v1.0/seurat/seurat4_no_removal_scrublet_cellranger6.R -i /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1/outs -o /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/seurat/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1 -s HT425B1-S1H1A3Y1N1Z1_1Bmn1_1 --scrublet /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/scrubblet/RNA --detailed
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/
Rscript /diskmnt/Projects/Users/austins2/tools/seurat4_no_removal_scrublet_sctv2_cellranger6.R -i /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1/outs/raw_feature_bc_matrix/ -o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1/ -s HT425B1-S1H1A3Y1N1Z1_1Bmn1_1 --scrublet /diskmnt/Projects/HTAN_BRCA_primary/Cellranger/scrubblet/RNA/ --detailed
# HT179C1-T1A3Y1N1
conda activate r4.2
Rscript /diskmnt/Projects/HTAN_analysis_2/ST_interactions/data_processing/seurat_multi/seurat_pipeline_v5.0.combo.R \
-s HT179C1-T1A3Y1N1 \
-d /diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/HT179C1-T1A3Y1N1 \
-m /diskmnt/Projects/Users/allakarpova/Tools/anaconda3/envs/signac/bin/macs2 \
-o /diskmnt/Projects/HTAN_analysis_2/TNP_CRC/seurat/combo_v5.0 \
--chrom_size /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt \
--prf_min 1000 --pct_min 15 --ns_max 5 --pc_first 2 --pc_num 50
# C3L-03372_CPT0275960013_2022-01-31 and C3N-00663_CPT0087730015_2021-03-08 were initially processed using the same code and pipeline published here:
# Liu, J. et al. Multi-scale signaling and tumor evolution in high-grade gliomas. Cancer Cell 42, 1217-1238.e19 (2024). https://doi.org/10.1016/j.ccell.2024.06.004 
# The above objects were each used for their own separate ongoing projects.
# These individual samples from each project relevant to this study were then subset and reprocessed in a uniform manner with the seurat 5.0.1 environment: `spatial_driver_co-detection_paper/r4.2_env_used_to_generate_probe_flanks_last_update_20240201.yml`
# The process for the uniform reporocessing of these samples is documented in the script `spatial_driver_co-detection_paper/snRNA_scRNA_preprocessing/Updating_snRNA_scRNA_to_seurat5.R`



