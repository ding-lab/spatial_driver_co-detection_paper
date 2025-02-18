#!/usr/bin/env bash
# this is a tsv file with lists of xenium variant targeting probe names: /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv
# activate the seurat 5 conda environment: conda activate seurat5
# HT227P1-S1H1L1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT227P1-S1H1L1U1__20240411__214815 \
-s HT227P1-S1H1L1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/
# HT242P1-S1H4L4U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \ \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT242P1-S1H4L4U1__20240411__214815 \
-s HT242P1-S1H4L4U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/
# HT243B1-S1H1A4U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT243B1-SH1A4U1__20240411__214815 \
-s HT243B1-S1H1A4U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/
# HT260C1-Th1K1L1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT260C1-Th1K1L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs \
-s HT260C1-Th1K1L1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT268B1-Th1H3L1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT268B1-Th1H3L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs \
-s HT268B1-Th1H3L1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0
# HT284P1-S1H1A1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__HT284P1-S1H1A1U1__20240417__211045 \
-s HT284P1-S1H1A1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# HT305B1-S1H5A1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__HT305B1-S1H5A1U1__20240501__171738 \
-s HT305B1-S1H5A1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# HT308B1-S1H5A4U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__HT308B1-S1H5A4U1__20240411__214816 \
-s HT308B1-S1H5A4U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# HT425B1-S1H1Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__HT425B1-S1H1Fp1U1__20240501__171738  \
-s HT425B1-S1H1Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# HT65B1-H1A1A4U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__HT65B1-H1A1A4U1__20240417__211045 \
-s HT65B1-H1A1A4U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE82256U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__WUPE82256U1-Fp1U1__20240411__214816 \
-s WUPE82256U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE25723U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE25723U1-Fp1U1__20240501__171738 \
-s WUPE25723U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE08092A1-S1H3U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__WUPE08092A1-S1H3U1__20240411__214816 \
-s WUPE08092A1-S1H3U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE68972U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__WUPE68972U1-Fp1U1__20240411__214816 \
-s WUPE68972U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE38758U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE38758U1-Fp1U1__20240501__171738 \
-s WUPE38758U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE62598U1-Fp2U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE62598U1-Fp2U1__20240501__171738 \
-s WUPE62598U1-Fp2U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE28077U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE28077U1-Fp1U1__20240501__171738 \
-s WUPE28077U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# WUPE04916U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__WUPE04916U1-Fp1U1__20240417__211045 \
-s WUPE04916U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# SP001P1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium_primary/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__SP001P1-Fp1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs \
-s SP001P1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0
# WUPE86238U1-Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__WUPE86238U1-Fp1U1__20240417__211045 \
-s WUPE86238U1-Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# SN105H1-Ma1Fp2-5U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__SN105H1-Ma1Fp2-5U1__20240501__171738 \
-s SN105H1-Ma1Fp2-5U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# SN106H1-Ma1Fp2-7U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__SN106H1-Ma1Fp2-7U1__20240501__171738 \
-s SN106H1-Ma1Fp2-7U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# C3L-01287-11Us2_1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium_primary/data/20240613__175032__20240613_Prostate_SNP/output-XETG00122__0033825__C3L-01287-11Us2_1__20240613__175113/ \
-s C3L-01287-11Us2_1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# C3L-03372-12Us3_1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium_primary/data/20240613__175032__20240613_Prostate_SNP/output-XETG00122__0033825__C3L-03372-12Us3_1__20240613__175113/ \
-s C3L-03372-12Us3_1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# C3N-00663-11U3
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium_primary/data/20240613__175032__20240613_Prostate_SNP/output-XETG00122__0033825__C3N-00663-11U3__20240613__175113/ \
-s C3N-00663-11U3 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects
# HT270P1-S1H1A1US2_1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__HT270P1-S1H1A1US2_1__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/  \
-s HT270P1-S1H1A1US2_1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0
# SP002C1-Fp1U2
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__SP002C1-Fp1U2__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s SP002C1-Fp1U2 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0
# HT060P1-S1R1Fp1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT060P1-S1R1Fp1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT060P1-S1R1Fp1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT061P1-S1P1A1L1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT061P1-S1P1A1L1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT061P1-S1P1A1L4U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L4U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT061P1-S1P1A1L4U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT125P1-S1H4A1L1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT125P1-S1H4A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT125P1-S1H4A1L1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT125P1-S1H8A1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT125P1-S1H8A1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT125P1-S1H8A1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT179C1-T1Fp3L5U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT179C1-T1Fp3L5U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT179C1-T1Fp3L5U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/
# HT185P1-S1H2L1U1
Rscript /diskmnt/Projects/Users/austins2/tools/seurat_5.0.1_xenium_v5.R \
-i /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT185P1-S1H2L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/ \
-s HT185P1-S1H2L1U1 \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_snv_panel_v1_snv_probe_names.tsv \
-o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/xr_v2_0/