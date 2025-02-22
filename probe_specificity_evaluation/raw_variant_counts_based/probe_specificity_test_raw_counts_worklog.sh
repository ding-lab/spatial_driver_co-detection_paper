# removed_unknown_v6/counts based
# checking genotyping results for Xenium SNV panel v1 samples (probe_specificity_test with probe based count excluding the unknown/necrotic cell types because they are more likely to have high amount of false positive binding)
# conda activate seurat5
# these rely on the Seurat objects as input.
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/
# 1
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/SP001P1-Fp1U1/SP001P1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/SP001P1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/SP001P1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s SP001P1-Fp1U1 \
-t 'Tumor' \
-u 'Necrosis'
# 2
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT242P1-S1H4L4U1/HT242P1-S1H4L4U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT242P1-S1H4L4U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT242P1-S1H4L4U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT242P1-S1H4L4U1 \
-t 'Tumor' \
-u 'LowCount'
# 3
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT227P1-S1H1L1U1/HT227P1-S1H1L1U1_processed_sub_2_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT227P1-S1H1L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT227P1-S1H1L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT227P1-S1H1L1U1 \
-t 'Tumor,PanIN' \
-u 'LowCount'
# 4
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT268B1-Th1H3L1U1/HT268B1-Th1H3L1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT268B1-Th1H3L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT268B1-Th1H3L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT268B1-Th1H3L1U1 \
-t 'Cancer cells proliferating,Cancer cells' \
-u 'Nothing'
# 5
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT243B1-S1H1A4U1/HT243B1-S1H1A4U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based//HT243B1-S1H1A4U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT243B1-S1H1A4U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT243B1-S1H1A4U1 \
-t 'Cancer cells,Cancer cells proliferating' \
-u 'Nothing'
# 6
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT308B1-S1H5A4U1/HT308B1-S1H5A4U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT308B1-S1H5A4U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT308B1-S1H5A4U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT308B1-S1H5A4U1 \
-t 'Cancer cells,Cancer cells proliferating' \
-u 'Nothing'
# 7
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT260C1-Th1K1L1U1/HT260C1-Th1K1L1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT260C1-Th1K1L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT260C1-Th1K1L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT260C1-Th1K1L1U1 \
-t 'Tumor_1,Tumor_2' \
-u 'Necrosis?,LowCount'
# 8
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE82256U1-Fp1U1/WUPE82256U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE82256U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE82256U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE82256U1-Fp1U1 \
-t 'PC' \
-u 'Unk_Low'
# 9
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE68972U1-Fp1U1/WUPE68972U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE68972U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE68972U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE68972U1-Fp1U1 \
-t 'PC,MKC_PC' \
-u 'Unk_CTSK,Unk_Low'
# 10
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT305B1-S1H5A1U1/HT305B1-S1H5A1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT305B1-S1H5A1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT305B1-S1H5A1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT305B1-S1H5A1U1  \
-t 'Cancer cells,Cancer cells proliferating' \
-u 'Nothing'
# 11
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT65B1-H1A1A4U1/HT65B1-H1A1A4U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT65B1-H1A1A4U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT65B1-H1A1A4U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT65B1-H1A1A4U1 \
-t 'Cancer cells,Cancer cells proliferating' \
-u ''
# 12
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE08092A1-S1H3U1/WUPE08092A1-S1H3U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE08092A1-S1H3U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE08092A1-S1H3U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE08092A1-S1H3U1 \
-t 'Tumor' \
-u 'Unknown'
# 13
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE62598U1-Fp2U1/WUPE62598U1-Fp2U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE62598U1-Fp2U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE62598U1-Fp2U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE62598U1-Fp2U1 \
-t 'Plasma2,Plasma1,B/Plasma2,Fibroblast/Adipo/Plasma_2,NK/T/Plasma1' \
-u ''
# 14
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/SN106H1-Ma1Fp2-7U1/SN106H1-Ma1Fp2-7U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/SN106H1-Ma1Fp2-7U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/SN106H1-Ma1Fp2-7U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s SN106H1-Ma1Fp2-7U1 \
-t 'B' \
-u 'Unknown_low_count'
# 15
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE28077U1-Fp1U1/WUPE28077U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE28077U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE28077U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE28077U1-Fp1U1 \
-t 'Plasma_2,B/Plasma2' \
-u ''
# 16
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE86238U1-Fp1U1/WUPE86238U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE86238U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE86238U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE86238U1-Fp1U1 \
-t 'B/Plasma_2/Plasma_1' \
-u ''
# 17
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT284P1-S1H1A1U1/HT284P1-S1H1A1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT284P1-S1H1A1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT284P1-S1H1A1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT284P1-S1H1A1U1 \
-t 'ITPN,Tumor_proliferating' \
-u 'Unknown,Necrosis?'
# 18
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT425B1-S1H1Fp1U1/HT425B1-S1H1Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT425B1-S1H1Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT425B1-S1H1Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT425B1-S1H1Fp1U1 \
-t 'Cancer cells' \
-u 'Nothing'
# 19
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE25723U1-Fp1U1/WUPE25723U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE25723U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE25723U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE25723U1-Fp1U1 \
-t 'Plasma2' \
-u 'Unknown'
# 20
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE38758U1-Fp1U1/WUPE38758U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE38758U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE38758U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE38758U1-Fp1U1 \
-t 'B/Plasma2,CCL5+/B/Plasma1' \
-u ''
# 21
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/WUPE04916U1-Fp1U1/WUPE04916U1-Fp1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/WUPE04916U1-Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/WUPE04916U1-Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s WUPE04916U1-Fp1U1 \
-t 'Fibroblast/Plasma2,Plasma2,B' \
-u 'Unknown'
# 22
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/SN105H1-Ma1Fp2-5U1/SN105H1-Ma1Fp2-5U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/SN105H1-Ma1Fp2-5U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/SN105H1-Ma1Fp2-5U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s SN105H1-Ma1Fp2-5U1 \
-t 'B/Plasma_2' \
-u 'Unknown_low_count'
# 23
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/C3L-01287-11Us2_1/C3L-01287-11Us2_1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/C3L-01287-11Us2_1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/C3L-01287-11Us2_1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s C3L-01287-11Us2_1 \
-t 'cancer cell,cancer cell region 2' \
-u ''
# 24
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/C3L-03372-12Us3_1/C3L-03372-12Us3_1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/C3L-03372-12Us3_1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/C3L-03372-12Us3_1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s C3L-03372-12Us3_1 \
-t 'Tumor,Tumor_prolif' \
-u 'LowCount'
# 25
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/C3N-00663-11U3/C3N-00663-11U3_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/C3N-00663-11U3/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/C3N-00663-11U3_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s C3N-00663-11U3 \
-t 'Tumor,Tumor_prolif' \
-u 'LowCount'
# 26
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT270P1-S1H1A1US2_1/HT270P1-S1H1A1US2_1_processed_sub_5_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT270P1-S1H1A1US2_1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT270P1-S1H1A1US2_1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_variant_allele_probes_hyphens.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_reference_allele_probes_hyphens.tsv \
-s HT270P1-S1H1A1US2_1 \
-t 'Tumor,PanIN' \
-u 'LowCount'
# 27
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/SP002C1-Fp1U2/SP002C1-Fp1U2_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/SP002C1-Fp1U2/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/SP002C1-Fp1U2_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_variant_allele_probes_hyphens.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_reference_allele_probes_hyphens.tsv \
-s SP002C1-Fp1U2 \
-t 'Tumor' \
-u 'Necrosis'
# 28
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT060P1-S1R1Fp1U1/HT060P1-S1R1Fp1U1_processed_sub.4.12.cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT060P1-S1R1Fp1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT060P1-S1R1Fp1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT060P1-S1R1Fp1U1 \
-t 'Tumor' \
-u 'LowCount'
# 29
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT061P1-S1P1A1L1U1/HT061P1-S1P1A1L1U1_processed_sub_4_6_8_12_14_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT061P1-S1P1A1L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT061P1-S1P1A1L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT061P1-S1P1A1L1U1 \
-t 'Tumor,PanIN' \
-u 'LowCount'
# 30
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT061P1-S1P1A1L4U1/HT061P1-S1P1A1L4U1_processed_sub_4_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT061P1-S1P1A1L4U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT061P1-S1P1A1L4U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT061P1-S1P1A1L4U1 \
-t 'Tumor' \
-u 'LowCount'
# 31
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT125P1-S1H4A1L1U1/HT125P1-S1H4A1L1U1_processed_sub_6_7_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT125P1-S1H4A1L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT125P1-S1H4A1L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT125P1-S1H4A1L1U1 \
-t 'Tumor' \
-u 'LowCount'
# 32
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT125P1-S1H8A1U1/HT125P1-S1H8A1U1_processed_sub_2_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT125P1-S1H8A1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT125P1-S1H8A1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT125P1-S1H8A1U1 \
-t 'Tumor' \
-u 'LowCount'
# 33
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT179C1-T1Fp3L5U1/HT179C1-T1Fp3L5U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT179C1-T1Fp3L5U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT179C1-T1Fp3L5U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT179C1-T1Fp3L5U1 \
-t 'Tumor' \
-u 'Unknown,Necrosis'
# 34
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT185P1-S1H2L1U1/HT185P1-S1H2L1U1_processed_sub_1_6_cluster.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT185P1-S1H2L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT185P1-S1H2L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT185P1-S1H2L1U1 \
-t 'Tumor' \
-u 'LowCount'

# now doing the removed_unknown_v6 counts based again but this time only including the cells thought to be tumor subclones as the cancer cell labels when using the `-t` flag.
# removed_unknown_v6_subclone/counts_based
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based
ln -s /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/* .
# for samples with subclones present remove the non-subclonal result links from the new folder since we want to preserve them and not overwrite it
rm C3L-01287-11Us2_1
rm HT260C1-Th1K1L1U1
rm probe_table_list.rds # scripts will generate a subclone-specific version of this output object 
rm Variant_specific_results # scripts will generate a subclone-specific version of this result folder 
rm variant_result_compiling_input_table_6.tsv # will need to make a subclone-specific version of this input table 
# 7 corrected run
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/HT260C1-Th1K1L1U1/HT260C1-Th1K1L1U1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/HT260C1-Th1K1L1U1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT260C1-Th1K1L1U1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s HT260C1-Th1K1L1U1 \
-t 'Tumor_2' \
-u 'Necrosis?,LowCount'
# 23 corrected run
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
-i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/C3L-01287-11Us2_1/C3L-01287-11Us2_1_processed.rds \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/C3L-01287-11Us2_1/ \
-c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/C3L-01287-11Us2_1_cell_types_v7.tsv \
-v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
-r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
-s C3L-01287-11Us2_1 \
-t 'cancer cell region 2' \
-u ''
# the output "Sample_ID_probe_specificity_table.csv" files from the above commands are listed in an input table called "variant_result_compiling_input_table_6.tsv" and used as input for the following steps.
# compiling individual table results to a single table that is just the results from individual samples.
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/
# the output of the following command is Supplementary Table 5: "Supp Table 5 count individual" 
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/variant_result_compiling_counts_v5.R variant_result_compiling_input_table_6.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/
# Supplementary Table 5 in a file called "All_variants_probe_specificity_results_by_sample.tsv" is used as input for the following command. The output of the following command is Supplementary Table 4: "Supp Table 4 count" and saved to the file called "All_variants_probe_specificity_results_removed_unknown_v6_subclone_counts_with_germline.tsv"
Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/probe_specificity_sample_summary_counts_germline_v5.R All_variants_probe_specificity_results_by_sample.tsv All_variants_probe_specificity_results_removed_unknown_v6_subclone_counts_with_germline.tsv &> runlog_probe_specificity_sample_summary_counts_germline_v6_subclone.txt #the output of this is Supplementary Table 4 