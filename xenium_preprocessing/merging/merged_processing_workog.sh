# worklog for the generation and cell typing of merged objects: primary PDAC Xenium and all snRNA.
# for plotting of figures related to merged analysis see the contents of the `spatial_driver_co-detection_paper/plotting_snippets/` folder
# single SCT - proceeding with the output of this object. - relied on umap.50PC
cd cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/
conda activate seurat5
# the following script relies on subset_opt function stored in the `spatial_driver_co-detection_paper/xenium_preprocessing/seurat_object_generation/subset_obj_seurat.R` file. 
# If you are trying to run this yourself you will need to update the hardcoded path to the file.
Rscript Merge_PDAC_primary_Xenium_20241204_single_layer.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/merge_input_table_20241204.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/ \
-p PDAC_merge_primary_KRAS --counts_cutoff 0 --pc_num_individual 30


# Updating cell types for primary PDAC cases based on the merged object clustering. 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/
conda activate seurat5
# this following script was then used to generate dotplots of cell type markers for the Seurat clusters in order to then go through and roughly annotate cell type/state of different clusters.
# For the PDAC samples the marker genes were from this file which contains the cell types and marker genes for PDAC from Supplementary Table 12: /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt
# It is a two column file of tab-separated values containing the "Cell type" and "Marker" columns Supplementary Table 12.
# This first few lines of this table look like the following:
# Gene_set	Gene
# Pancreatic_epithelial	EHF
# Pancreatic_epithelial	EPCAM
# Acinar	AMY2A
# Acinar	GATM
# Acinar	AQP8
# Duct_like_1	CFTR
# Duct_like_1	FXYD2
# by changing out the markerset for others listed in Supplementary Table 12 we used a similar method for plotting markers in other individual objects.
# the script run in the following commands is located at `spatial_driver_co-detection_paper/xenium_preprocessing/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R`
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_seurat_cluster.csv \
seurat_clusters
# Following an initial first pass cell typing, clusters were then subclustered. The subclusters were then annotated based on detection of marker geneset.
# The subclusters were then annotated based on their expression of the same marker genesets.
# the record of how subclusters were calculated is recorded in `spatial_driver_co-detection_paper/merging/xenium/PDAC_merge_primary_KRAS_subclustering.R`
# While inside of this R session documented in `PDAC_merge_primary_KRAS_subclustering.R` I also calculated the degs for different subcluster results using and annotated which ones were shared across what clusters using system() calls to the script here
# `spatial_driver_co-detection_paper/xenium_preprocessing/Shared_marker_comparison_from_seurat_degs.R` . This script specifically works with the output of the FindAllMarkers() function from Seurat.
# During the subclustering of each cluster/group of clusters in the Xenium object the above script was then used again to generate a different set of dotplots based on the subclustered labels.
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/
conda activate seurat5
# cell typing subclusters under the merge.sub.19 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_merge.sub.19.tsv \
merge.sub.19
# cell typing subclusters under the merge.sub.13 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_merge.sub.13.tsv \
merge.sub.13
# cell typing subclusters under the unified.merge.sub.2.20.19_1.19_5 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.2.20.19_1.19_5.tsv \
unified.merge.sub.2.20.19_1.19_5
# cell typing subclusters under the unified.merge.sub.1.14.25.19_0.19_4.19_6 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.1.14.25.19_0.19_4.19_6.tsv \
unified.merge.sub.1.14.25.19_0.19_4.19_6
# cell typing subclusters under the unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.tsv \
unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8
# cell typing subclusters under the unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4.tsv \
unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4
# cell typing subclusters under the unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_5 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_5.tsv \
unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_5
# cell typing subclusters under the unified.merge.sub.0.21.19_2.19_7.13_3.13_0 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.tsv \
unified.merge.sub.0.21.19_2.19_7.13_3.13_0
# cell typing subclusters under the unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4.tsv \
unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4
# cell typing subclusters under the unified.merge.sub.3.11.res_0.3 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.3.11.res_0.3.tsv \
unified.merge.sub.3.11.res_0.3
# cell typing subclusters under the unified.merge.sub.6.10.res_0.3 column in the metadata
Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.6.10.res_0.3.tsv \
unified.merge.sub.6.10.res_0.3

# making the merged snRNA object
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/sample_mutation_table_based_on_bulk_7_with_snRNA.tsv
# conda activate seurat5
Rscript Merge_all_samples_snRNA_20250119_separate_layers.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/merge_input_table_snRNA_20250119.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/ \
-p All_samples_merge_snRNA_20250119_separateSCT --counts_cutoff 0 --pc_num_individual 30
# Retained the most up-to-date versions of cell type labels from prior work for each of the snRNA objects. These were included in the individual snRNA ojects.
# An example of this can be found in the `spatial_driver_co-detection_paper/plotting_snippets/plotting_snippets.R` file for how this pertained to HT065B1

