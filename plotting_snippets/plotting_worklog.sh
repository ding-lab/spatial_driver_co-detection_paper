# plotting all Xenium object genotypes
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/
conda activate seurat5
# These featureplot outputs of the below command were evaluated along with the outputs of the `spatial_driver_co-detection_paper/probe_specificity_evaluation/`, `spatial_driver_co-detection_paper/plotting_snippets/vaf_barplots_cancer_vs_normal/`, `spatial_driver_co-detection_paper/plotting_snippets/vaf_barplots_cancer_vs_normal/`
# This command uses the script and input table from `spatial_driver_co-detection_paper/spatial_and_umap_with_features_and_label/plotting_all_variant_probes/` . It also relies on auxillary files from `spatial_driver_co-detection_paper/probe_specificity_evaluation/`
# The featureplots output from this command were used as part of the manual review process and Extended Data Figure 2b
Rscript /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/plotting_genotypes.R \
plotting_genotypes_input_table.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/

# the file `spatial_driver_co-detection_paper/neo_norm_unk_table_v7.tsv` was an auxillary table with paths to object files and cell type labels which was used as an input table for the following command includes several columns:
# # Sample_ID: the sample ID
# # Xenium_neoplastic_labels_v6: the labels of all neoplastic cells in the sample
# # Xenium_neoplastic_labels_v6_subclone: the labels of all of the neoplastic cells in the sample that belong to the clone with the mutation. If the mutation has not been found to be subclonal then this column is identical to Xenium_neoplastic_labels_v6
# # Xenium_unknown_low_quality_labels_v6: These are the labels of the low-quality and/or low-count cells that could not be accurately determined to be neoplastic cells or normal cells based on their morphology and expression profiles. These cells were subset out and excluded during the manual review process.
# # Manual_xenium_cell_types_v6: These are tsv files containing the cell types of the sample. The first column is the cell barcode. The second column is the cell type/state.
# # Xenium_snv_object: This is the file for the seurat object.
# # feature_csv: This is a list of all variant site probes (variant and reference) for which the sample was found to be mutated following manual review. In the case of a 3-way probe design all all alleles for the targeted site are listed.
# # mutant_csv: this is the same as the feature csv but only contains the probes targeting the alleles present in the sample in the case of the 3-way probe design.
# # feature_csv_table_plots: duplicate column of feature_csv column.

# plotting the updated neoplastic normal unknown labels (output file pattern is *_neoplastic_clone_normal_unknown_v4*pdf or *_neoplastic_normal_unknown_v4*pdf) for the Xenium spatial driver project
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/
conda activate seurat5
mkdir darkblue_ceiling1_v7
# The following 2 scripts rely upon the following input tables:
# # Supplementary Table 5: "Supp Table 5 count individual" referred to as "all_sample_summary" in the script `spatial_driver_co-detection_paper/vaf_barplots_cancer_vs_normal/plotting_barplots_alt_count_v6.R`
# # Supplementary Table 8: "Supp Table 8 cell individual" referred to as "all_sample_summary_cell" in the script `spatial_driver_co-detection_paper/vaf_barplots_cancer_vs_normal/plotting_barplots_alt_count_v6.R` and "all_sample_summary" in the script `spatial_driver_co-detection_paper/vaf_barplots_cancer_vs_normal/plotting_barplots_alt_count_v6.R`
# # the `spatial_driver_co-detection_paper/neo_norm_unk_table_v7.tsv` table referred to as `input_table` in both `spatial_driver_co-detection_paper/vaf_barplots_cancer_vs_normal/plotting_barplots_alt_count_v6.R` and `spatial_driver_co-detection_paper/vaf_barplots_cancer_vs_normal/plotting_barplots_alt_count_v6.R`
Rscript plotting_barplots_alt_cells_v6.R &> plotting_barplots_alt_cells_v6.log
Rscript plotting_barplots_alt_count_v6.R &> plotting_barplots_alt_count_v6.log
# The following script relies upon the these input tables:
# # Supplementary Table 8: "Supp Table 8 cell individual" referred to as "all_sample_summary_cell" in the script `spatial_driver_co-detection_paper/vaf_boxplots_cancer_vs_normal/plotting_boxplots_alt_cells_false_positive_v6_subclone.R`
# # the `spatial_driver_co-detection_paper/neo_norm_unk_table_v7.tsv` table referred to as `input_table` in `spatial_driver_co-detection_paper/vaf_boxplots_cancer_vs_normal/plotting_boxplots_alt_cells_false_positive_v6_subclone.R`
Rscript plotting_boxplots_alt_cells_false_positive_v6_subclone.R # This generates one PDF per sample in the `input_table` and eahc PDF contains a box plot for every variant probe present in the value of the `feature_csv` column for that sample. The output plots are included on Extended Data Figure 3b and 3c, Extended Data Figure 6e.

# the `spatial_driver_co-detection_paper/neo_norm_unk_table_v7.tsv` table used as input for the following command includes several columns:
# # Sample_ID: the sample ID
# # Xenium_neoplastic_labels_v6: the labels of all neoplastic cells in the sample
# # Xenium_neoplastic_labels_v6_subclone: the labels of all of the neoplastic cells in the sample that belong to the clone with the mutation. If the mutation has not been found to be subclonal then this column is identical to Xenium_neoplastic_labels_v6
# # Xenium_unknown_low_quality_labels_v6: These are the labels of the low-quality and/or low-count cells that could not be accurately determined to be neoplastic cells or normal cells based on their morphology and expression profiles. These cells were subset out and excluded during the manual review process.
# # Manual_xenium_cell_types_v6: These are tsv files containing the cell types of the sample. The first column is the cell barcode. The second column is the cell type/state.
# # Xenium_snv_object: This is the file for the seurat object.
# # feature_csv: This is a list of all variant site probes (variant and reference) for which the sample was found to be mutated following manual review. In the case of a 3-way probe design all all alleles for the targeted site are listed.
# # mutant_csv: this is the same as the feature csv but only contains the probes targeting the alleles present in the sample in the case of the 3-way probe design.
# # feature_csv_table_plots: duplicate column of feature_csv column.
Rscript plotting_neoplastic_labels_7_darkblue_ceiling1.R neo_norm_unk_table_v7.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/
# this script generates output plots for each sample in the input table. Plots from this script were used in figure 3e, Figure 4c, Extended Data figure 2a, and Extended Data figure 6d.

conda activate seurat5
plotting_snippets.R # this file was run interactively in an R session and records how various plots using in figure generation were generated. 