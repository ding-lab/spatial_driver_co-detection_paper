#!/usr/bin/env bash
# commands to resegment xenium datasets with uniform cell segmentation parameters across all datasets
# Add Xenium ranger to your path in ~/.bashrc file before executing: export PATH=/diskmnt/Projects/Users/austins2/software/xeniumranger-xenium2.0.0:$PATH
cd /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/
# HT260C1-Th1K1L1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT260C1-Th1K1L1U1__20240322__175654 --id output-XETG00122__0010297__HT260C1-Th1K1L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true

# HT268B1-Th1H3L1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT268B1-Th1H3L1U1__20240322__175654 --id output-XETG00122__0010297__HT268B1-Th1H3L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true

# SP001P1-Fp1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__SP001P1-Fp1U1__20240322__175654 --id output-XETG00122__0010297__SP001P1-Fp1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true

cd /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/
# HT270P1-S1H1A1US2_1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__HT270P1-S1H1A1US2_1__20240112__230102 --id output-XETG00122__0010375__HT270P1-S1H1A1US2_1__20240112__230102_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true

# SP002C1-Fp1U2
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__SP002C1-Fp1U2__20240112__230102 --id output-XETG00122__0010375__SP002C1-Fp1U2__20240112__230102_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true

cd /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/
# HT060P1-S1R1Fp1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT060P1-S1R1Fp1U1__20241108__215027 --id output-XETG00122__0033854__HT060P1-S1R1Fp1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable

# HT061P1-S1P1A1L1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L1U1__20241108__215027 --id output-XETG00122__0034249__HT061P1-S1P1A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable

# HT061P1-S1P1A1L4U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L4U1__20241108__215027 --id output-XETG00122__0034249__HT061P1-S1P1A1L4U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable

# HT125P1-S1H4A1L1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT125P1-S1H4A1L1U1__20241108__215027 --id output-XETG00122__0033854__HT125P1-S1H4A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable

# HT125P1-S1H8A1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT125P1-S1H8A1U1__20241108__215027 --id output-XETG00122__0034249__HT125P1-S1H8A1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable

# HT179C1-T1Fp3L5U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT179C1-T1Fp3L5U1__20241108__215027 --id output-XETG00122__0033854__HT179C1-T1Fp3L5U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable

# HT185P1-S1H2L1U1
xeniumranger resegment --xenium-bundle /diskmnt/primary/Xenium/data/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT185P1-S1H2L1U1__20241108__215027 --id output-XETG00122__0034249__HT185P1-S1H2L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true --expansion-distance 5 --dapi-filter 100 --resegment-nuclei true --localcores 20 --interior-stain disable --boundary-stain disable
