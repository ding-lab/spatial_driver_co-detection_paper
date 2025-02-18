#!/usr/bin/env bash
# This is the worklog for alignment of H&E images generated on post xenium slides to the Xenium dapi image 
# Each image was aligned using at least 3 settings. 
# The resulting he_aligned.ome.tif image that was generated from each run was then loaded into Xenium Explorer to qualitatively evaluate which of the aligned images was optimal for interpretation.
# This was done by looking at how well the Xenium nuclei borders in Xenium Explorer overlapped the nuclei in the H&E in the center of the tissue and at each of the corners of the image.
# activate the python 3.10 conda environment from HEX-SIFT: conda activate py3.10
# some notes: 
#   - there is no difference in the output between image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py and image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py  The only difference is that image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py can read the images written by the Xenium Onboard Analysis software v3.0+ as input while image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py cannot.
# HT227P1-S1H1L1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT227P1-S1H1L1U1
ln -s /diskmnt/primary/Xenium_primary/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT227P1-S1H1L1U1__20240411__214815/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT227P1-S1H1L1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT227P1-S1H1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT227P1-S1H1L1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT227P1-S1H1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT227P1-S1H1L1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT227P1-S1H1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT242P1-S1H4L4U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT242P1-S1H4L4U1
ln -s /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT242P1-S1H4L4U1__20240411__214815/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT242P1-S1H4L4U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT242P1-S1H4L4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT242P1-S1H4L4U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT242P1-S1H4L4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT242P1-S1H4L4U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT242P1-S1H4L4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT243B1-S1H1A4U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT243B1-S1H1A4U1
ln -s /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022792__HT243B1-SH1A4U1__20240411__214815/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT243B1-S1H1A4U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT243B1-S1H1A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT243B1-S1H1A4U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT243B1-S1H1A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT243B1-S1H1A4U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT243B1-S1H1A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT260C1-Th1K1L1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT260C1-Th1K1L1U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT260C1-Th1K1L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT260C1-Th1K1L1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT260C1-Th1K1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT260C1-Th1K1L1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT260C1-Th1K1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT260C1-Th1K1L1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT260C1-Th1K1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT268B1-Th1H3L1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT268B1-Th1H3L1U1
ln -s /diskmnt/primary/Xenium_primary/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__HT268B1-Th1H3L1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT268B1-Th1H3L1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT268B1-Th1H3L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT268B1-Th1H3L1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT268B1-Th1H3L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT268B1-Th1H3L1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT268B1-Th1H3L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
mkdir v9_bs0m100
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT268B1-Th1H3L1U1/v9_bs0m100
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT268B1-Th1H3L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 100
# HT284P1-S1H1A1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT284P1-S1H1A1U1
ln -s /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__HT284P1-S1H1A1U1__20240417__211045/morphology_focus/morphology_focus_0000.ome.tif .
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT284P1-S1H1A1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT284P1-S1H1A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT284P1-S1H1A1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT284P1-S1H1A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT284P1-S1H1A1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT284P1-S1H1A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT305B1-S1H5A1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT305B1-S1H5A1U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__HT305B1-S1H5A1U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT305B1-S1H5A1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT305B1-S1H5A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT305B1-S1H5A1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT305B1-S1H5A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT305B1-S1H5A1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT305B1-S1H5A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT308B1-S1H5A4U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT308B1-S1H5A4U1
ln -s /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__HT308B1-S1H5A4U1__20240411__214816/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT308B1-S1H5A4U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT308B1-S1H5A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT308B1-S1H5A4U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT308B1-S1H5A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT308B1-S1H5A4U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT308B1-S1H5A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT425B1-S1H1Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT425B1-S1H1Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__HT425B1-S1H1Fp1U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT425B1-S1H1Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT425B1-S1H1Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT425B1-S1H1Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT425B1-S1H1Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT425B1-S1H1Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT425B1-S1H1Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT65B1-H1A1A4U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT65B1-H1A1A4U1
ln -s /diskmnt/primary/Xenium_primary/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__HT65B1-H1A1A4U1__20240417__211045/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT65B1-H1A1A4U1/v9_bs0m60 
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT65B1-H1A1A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT65B1-H1A1A4U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT65B1-H1A1A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/HT65B1-H1A1A4U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../HT65B1-H1A1A4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# WUPE82256U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE82256U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__WUPE82256U1-Fp1U1__20240411__214816/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE82256U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE82256U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE82256U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE82256U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE82256U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE82256U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
mdkir v9_bs0m100
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE82256U1-Fp1U1/v9_bs0m100
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE82256U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100
# WUPE25723U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE25723U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE25723U1-Fp1U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE25723U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE25723U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE25723U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE25723U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE25723U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE25723U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# WUPE08092A1-S1H3U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE08092A1-S1H3U1
ln -s /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__WUPE08092A1-S1H3U1__20240411__214816/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE08092A1-S1H3U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE08092A1-S1H3U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE08092A1-S1H3U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE08092A1-S1H3U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE08092A1-S1H3U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE08092A1-S1H3U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# WUPE68972U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE68972U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240411__214732__20240411_PECGS_HTAN_SNP-eval/output-XETG00122__0022898__WUPE68972U1-Fp1U1__20240411__214816/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE68972U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE68972U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE68972U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE68972U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE68972U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE68972U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20 
# WUPE38758U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE38758U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE38758U1-Fp1U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE38758U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE38758U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE38758U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE38758U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE38758U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE38758U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20 
# WUPE62598U1-Fp2U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE62598U1-Fp2U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE62598U1-Fp2U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE62598U1-Fp2U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE62598U1-Fp2U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE62598U1-Fp2U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE62598U1-Fp2U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE62598U1-Fp2U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE62598U1-Fp2U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# WUPE28077U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE28077U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022793__WUPE28077U1-Fp1U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE28077U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE28077U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE28077U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE28077U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE28077U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE28077U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# WUPE04916U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE04916U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__WUPE04916U1-Fp1U1__20240417__211045/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE04916U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE04916U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE04916U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE04916U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE04916U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE04916U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20 # this one was the best and used in subsequent manual review/analysis/figures
# SP001P1-Fp1U1
# aligning image for SP001P1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP001P1-Fp1U1/
ln -s /diskmnt/primary/Xenium_primary/data/data_reprocess_xeniumranger_v2.0/20240322__175618__202040322_HTAN_SNP_prostate/output-XETG00122__0010297__SP001P1-Fp1U1__20240322__175654_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_0000.ome.tif
mdkir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP001P1-Fp1U1/v9_bs4m20
python3 $austin/tools/xenium_image_alignment/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SP001P1-Fp1U1_input.ome.tif -d ../morphology_mip.ome.tif -b -s 4 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP001P1-Fp1U1/v9_bs2m40
python3 $austin/tools/xenium_image_alignment/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SP001P1-Fp1U1_input.ome.tif -d ../morphology_mip.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP001P1-Fp1U1/v9_bs0m60
python3 $austin/tools/xenium_image_alignment/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SP001P1-Fp1U1_input.ome.tif -d ../morphology_mip.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
# WUPE86238U1-Fp1U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE86238U1-Fp1U1
ln -s /diskmnt/primary/Xenium/data/20240417__211003__20240417_IRD_SNP/output-XETG00122__0033835__WUPE86238U1-Fp1U1__20240417__211045/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE86238U1-Fp1U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE86238U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -f -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE86238U1-Fp1U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE86238U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -f -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/WUPE86238U1-Fp1U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../WUPE86238U1-Fp1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -f -s 4 -m 20
# SN105H1-Ma1Fp2-5U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN105H1-Ma1Fp2-5U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__SN105H1-Ma1Fp2-5U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN105H1-Ma1Fp2-5U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SN105H1-Ma1Fp2-5U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN105H1-Ma1Fp2-5U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SN105H1-Ma1Fp2-5U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN105H1-Ma1Fp2-5U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SN105H1-Ma1Fp2-5U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# SN106H1-Ma1Fp2-7U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN106H1-Ma1Fp2-7U1
ln -s /diskmnt/primary/Xenium/data/20240501__171638__20240501_SenNet_HTAN_PECGS_SNP/output-XETG00122__0022794__SN106H1-Ma1Fp2-7U1__20240501__171738/morphology_focus/morphology_focus_0000.ome.tif
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN106H1-Ma1Fp2-7U1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SN106H1-Ma1Fp2-7U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN106H1-Ma1Fp2-7U1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SN106H1-Ma1Fp2-7U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SN106H1-Ma1Fp2-7U1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_pre_xoa_v3.0.py -i ../SN106H1-Ma1Fp2-7U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# I had to update HEX-SIFT to read the input files from samples that were stained with the multimodal protein stain. 
# The underlying keypoint/descriptor calculation, matching calculation, and perspective warps are the same. you can check this in my commit history on the HEX-SIFT github.
# Samples with multimodal protein stains had their H&E images aligned with image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py because it can handle the new input image format.
# C3L-01287-11Us2_1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1
ln -s /diskmnt/primary/Xenium_primary/data/20240613__175032__20240613_Prostate_SNP/output-XETG00122__0033825__C3L-01287-11Us2_1__20240613__175113/morphology_focus/morphology_focus_000* .
mkdir v10_bs0m60 v10_bs0m60_multimodal v10_bs2m40 v10_bs4m20 v10_bs4m20_multimodal v10_bs4m15_multimodal v10_bs5m20 v10_bs5m20_multimodal v10_bs5m10 v10_s4m20_multimodal v10_s0m60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs0m60 
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs0m60_multimodal
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif --multimodal_path ../ -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs4m20_multimodal
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif --multimodal_path ../ -b -s 4 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs4m15_multimodal
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif --multimodal_path ../ -b -s 4 -m 15
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs5m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 5 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs5m20_multimodal
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif --multimodal_path ../ -b -s 5 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_bs5m10
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 5 -m 15
# non-blue
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_s4m20_multimodal
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif --multimodal_path ../ -s 4 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-01287-11Us2_1/v10_s0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-01287-11Us2_1.ome.tif --multimodal_path ../ -s 0 -m 60
# C3L-03372-12Us3_1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-03372-12Us3_1
ln -s /diskmnt/primary/Xenium_primary/data/20240613__175032__20240613_Prostate_SNP/output-XETG00122__0033825__C3L-03372-12Us3_1__20240613__175113/morphology_focus/morphology_focus_000* .
mkdir v9_bs0m60 v9_bs2m40 v9_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-03372-12Us3_1/v9_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-03372-12Us3_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-03372-12Us3_1/v9_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-03372-12Us3_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3L-03372-12Us3_1/v9_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3L-03372-12Us3_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# C3N-00663-11U3
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3N-00663-11U3
ln -s /diskmnt/primary/Xenium_primary/data/20240613__175032__20240613_Prostate_SNP/output-XETG00122__0033825__C3N-00663-11U3__20240613__175113/morphology_focus/morphology_focus_000* .
mkdir v10_bs0m60 v10_bs2m40 v10_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3N-00663-11U3/v10_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3N-00663-11U3.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3N-00663-11U3/v10_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3N-00663-11U3.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/C3N-00663-11U3/v10_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../C3N-00663-11U3.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20 
# HT270P1-S1H1A1US2_1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT270P1-S1H1A1US2_1
ln -s /diskmnt/primary/Xenium_primary/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__HT270P1-S1H1A1US2_1__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_0000.ome.tif
mkdir v10_bs0m100 v10_bs0m20 v10_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT270P1-S1H1A1US2_1/v10_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT270P1-S1H1A1US2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT270P1-S1H1A1US2_1/v10_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT270P1-S1H1A1US2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT270P1-S1H1A1US2_1/v10_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT270P1-S1H1A1US2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT270P1-S1H1A1US2_1/v10_bs0m100
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT270P1-S1H1A1US2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT270P1-S1H1A1US2_1/v10_bs0m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT270P1-S1H1A1US2_1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 20
# SP002C1-Fp1U2
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/non-pdac/SP002C1-Fp1U2 
ln -s /diskmnt/primary/Xenium_primary/data/data_reprocess_xeniumranger_v2.0/20240112__230030__20240112_24018_8UTB89_retry/output-XETG00122__0010375__SP002C1-Fp1U2__20240112__230102_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_0000.ome.tif
mkdir v10_bs0m60s v10_bs2m40 v10_bs4m20 v10_bs0m100 v10_bs0m150
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP002C1-Fp1U2_v10_alignment/v10_bs0m60s
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../SP002C1-Fp1U2.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP002C1-Fp1U2_v10_alignment/v10_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../SP002C1-Fp1U2.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP002C1-Fp1U2_v10_alignment/v10_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../SP002C1-Fp1U2.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP002C1-Fp1U2_v10_alignment/v10_bs0m100
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../SP002C1-Fp1U2.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/SP002C1-Fp1U2_v10_alignment/v10_bs0m150
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../SP002C1-Fp1U2.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 150
# HT060P1-S1R1Fp1U1
# due to the placement of this tissue on the slide it was impossible to crop the scan of the H&E slide to just this piece of tissue if restricted to a square box so I wrote a script that let me crop an image with an irregularly shaped border
conda activate py3.10
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT060P1-S1R1Fp1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
# this command makes an all black image that is the same size as the input image specified with -i
python3 $austin/tools/HEX-SIFT/write_annotation_masked_image_v2.py -b -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/testing_blank_crop/ID_0033854_Scan1.qptiff -n HT060P1-S1R1Fp1U1_blank_image_for_annotation -d /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/
# I then loaded the original H&E slide scan into QuPath and drew an annotation around the entire piece of tissue for HT060P1-S1R1Fp1U1 from slide ID_0033854_Scan1.qptiff that did not contain tissue from any other samples.
# Load the all black image and the irregular shape annotation (qpdata file) (need to change the annotation color to white) of the region containing the tissue. 
# Export the image from QuPath using export > rendered RGB option (rendered annotation saved as file: "HT060P1-S1R1Fp1U1_blank_annotated.ome.tif") This generates an all black image with the region that will be cropped out of the actual image outlined in white. upload that to the cluster.
python3 $austin/tools/HEX-SIFT/write_annotation_masked_image_v2.py -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/testing_blank_crop/ID_0033854_Scan1.qptiff -a HT060P1-S1R1Fp1U1_blank_annotated.ome.tif -n HT060P1-S1R1Fp1U1_crop -d /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/ -w 10
# the HT060P1-S1R1Fp1U1_crop_white_background.ome.tif file that is generated from the above command was then loaded into QuPath and exported as an OME TIF with the following settings:
# compression type: ZLIB (lossless)
# pyramidal downsample: 2.0
# Tile size: 1024 px
# This is necessary because the image that is written by the above script is a little broken in how the layers of the image are saved. 
# HEX-SIFT expects the channels and dimensions of the image to be saved in a very specific format/order and when it writes the file it uses a sligthly different format so that it can write only a single file as output that is compatible with XeniumExplorer instead of 3. 
# I used the format I did for HEX-SIFT because this was the only format that would let me load the output aligned images into QuPath, Xenium Explorer, and ImageJ/Fiji, which makes working with pathologist annotations and repeatability of evaluating anything on the aligned images a lot simpler.
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/post_xoa_v3_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT060P1-S1R1Fp1U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/post_xoa_v3_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT060P1-S1R1Fp1U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/post_xoa_v3_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT060P1-S1R1Fp1U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
mkdir /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/post_xoa_v3_bs0m100
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT060P1-S1R1Fp1U1/post_xoa_v3_bs0m100
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT060P1-S1R1Fp1U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100 # this one was the best and used in subsequent manual review/analysis/figures
# HT061P1-S1P1A1L1U1
conda activate py3.10
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
# this command makes an all black image that is the same size as the input image specified with -i
python3 $austin/tools/HEX-SIFT/write_annotation_masked_image_v2.py -b -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/testing_blank_crop/ID_0034249_Scan1.qptiff -n HT061P1-S1P1A1L1U1_blank_image_for_annotation -d /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/
# I then loaded the original H&E slide scan into QuPath and drew an annotation around the entire piece of tissue for HT061P1-S1P1A1L1U1 from slide ID_0034249_Scan1.qptiff that did not contain tissue from any other samples.
# Load the all black image and the irregular shape annotation (qpdata file) (need to change the annotation color to white) of the region containing the tissue. 
# Export the image from QuPath using export > rendered RGB option (rendered annotation saved as file: "HT061P1-S1P1A1L1U1_blank_annotated_2.ome.tif") This generates an all black image with the region that will be cropped out of the actual image outlined in white. upload that to the cluster.
python3 $austin/tools/HEX-SIFT/write_annotation_masked_image_v2.py -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/testing_blank_crop/ID_0034249_Scan1.qptiff -a HT061P1-S1P1A1L1U1_blank_annotated_2.ome.tif -n HT061P1-S1P1A1L1U1_crop_2 -d /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/ -w 10
# the HT061P1-S1P1A1L1U1_crop_2_white_background.ome.tif file that is generated from the above command was then loaded into QuPath and exported as an OME TIF with the following settings:
# compression type: ZLIB (lossless)
# pyramidal downsample: 2.0
# Tile size: 1024 px
# This is necessary because the image that is written by the above script is a little broken in how the layers of the image are saved. 
# HEX-SIFT expects the channels and dimensions of the image to be saved in a very specific format/order and when it writes the file it uses a sligthly different format so that it can write only a single file as output that is compatible with XeniumExplorer instead of 3. 
# I used the format I did for HEX-SIFT because this was the only format that would let me load the output aligned images into QuPath, Xenium Explorer, and ImageJ/Fiji, which makes working with pathologist annotations and repeatability of evaluating anything on the aligned images a lot simpler.
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/post_xoa_v3_bs0m60 
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/post_xoa_v3_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/post_xoa_v3_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
mkdir /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/post_xoa_v3_bs0m100
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/post_xoa_v3_bs0m100 # used this one. did not end up needing the cropped image 
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100 # this one was the best and used in subsequent manual review/analysis/figures
mkdir /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/crop_2_post_xoa_v3_bs0m100
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L1U1/crop_2_post_xoa_v3_bs0m100
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L1U1_crop_2_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100
# HT061P1-S1P1A1L4U1
conda activate py3.10
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L4U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT061P1-S1P1A1L4U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L4U1/post_xoa_v3_bs0m60
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L4U1/post_xoa_v3_bs2m40
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT061P1-S1P1A1L4U1/post_xoa_v3_bs4m20
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT061P1-S1P1A1L4U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20 
# HT125P1-S1H4A1L1U1
conda activate py3.10
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H4A1L1U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT125P1-S1H4A1L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H4A1L1U1/post_xoa_v3_bs0m60 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT125P1-S1H4A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H4A1L1U1/post_xoa_v3_bs2m40 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT125P1-S1H4A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H4A1L1U1/post_xoa_v3_bs4m20 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT125P1-S1H4A1L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT125P1-S1H8A1U1
conda activate py3.10
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H8A1U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT125P1-S1H8A1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H8A1U1/post_xoa_v3_bs0m60 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT125P1-S1H8A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H8A1U1/post_xoa_v3_bs2m40 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT125P1-S1H8A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT125P1-S1H8A1U1/post_xoa_v3_bs4m20 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT125P1-S1H8A1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
# HT179C1-T1Fp3L5U1
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0033854__HT179C1-T1Fp3L5U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
# this command makes an all black image that is the same size as the input image specified with -i
python3 $austin/tools/HEX-SIFT/write_annotation_masked_image_v2.py -b -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/testing_blank_crop/ID_0033854_Scan1.qptiff -n HT179C1-T1Fp3L5U1_blank_image_for_annotation -d /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/
# I then loaded the original H&E slide scan into QuPath and drew an annotation around the entire piece of tissue for HT061P1-S1P1A1L1U1 from slide ID_0034249_Scan1.qptiff that did not contain tissue from any other samples.
# Load the all black image and the irregular shape annotation (qpdata file) (need to change the annotation color to white) of the region containing the tissue. 
# Export the image from QuPath using export > rendered RGB option (rendered annotation saved as file: "HT061P1-S1P1A1L1U1_blank_annotated_2.ome.tif") This generates an all black image with the region that will be cropped out of the actual image outlined in white. upload that to the cluster.
python3 $austin/tools/HEX-SIFT/write_annotation_masked_image_v2.py -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/testing_blank_crop/ID_0033854_Scan1.qptiff -a HT179C1-T1Fp3L5U1_blank_annotated.ome.tif -n HT179C1-T1Fp3L5U1_crop -d /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/ -w 10
# the HT179C1-T1Fp3L5U1_crop_white_background.ome.tif file that is generated from the above command was then loaded into QuPath and exported as an OME TIF with the following settings:
# compression type: ZLIB (lossless)
# pyramidal downsample: 2.0
# Tile size: 1024 px
# This is necessary because the image that is written by the above script "write_annotation_masked_image_v2.py" is a little broken in how the order of the layers of the image are saved. 
# HEX-SIFT expects the channels and dimensions of the image to be saved in a very specific format/order and when it writes the file it uses a sligthly different format so that it can write only a single file as output that is compatible with XeniumExplorer instead of 3. 
# I used the format I did for HEX-SIFT because this was the only format that would let me load the output aligned images into QuPath, Xenium Explorer, and ImageJ/Fiji, which makes working with pathologist annotations and repeatability of evaluating anything on the aligned images a lot simpler.
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/post_xoa_v3_bs0m60 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT179C1-T1Fp3L5U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/post_xoa_v3_bs2m40 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT179C1-T1Fp3L5U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/post_xoa_v3_bs4m20 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT179C1-T1Fp3L5U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20
mkdir /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/post_xoa_v3_bs0m100
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT179C1-T1Fp3L5U1/post_xoa_v3_bs0m100 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT179C1-T1Fp3L5U1_crop_white_background_fix.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 100
# HT185P1-S1H2L1U1
conda activate py3.10
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT185P1-S1H2L1U1
ln -s /diskmnt/primary/Xenium/data/data_reprocess_xeniumranger_v2.0/20241108__214933__20241108_V1_KRAS_HTAN_PDAC_CRC/output-XETG00122__0034249__HT185P1-S1H2L1U1__20241108__215027_xrv2_0_ne_5um_df_100_rn_true/outs/morphology_focus/morphology_focus_000* .
mkdir post_xoa_v3_bs0m60 post_xoa_v3_bs2m40 post_xoa_v3_bs4m20
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT185P1-S1H2L1U1/post_xoa_v3_bs0m60 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT185P1-S1H2L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 0 -m 60 # this one was the best and used in subsequent manual review/analysis/figures
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT185P1-S1H2L1U1/post_xoa_v3_bs2m40 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT185P1-S1H2L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 2 -m 40
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/image_alignment/HT185P1-S1H2L1U1/post_xoa_v3_bs4m20 #run
python3 $austin/tools/HEX-SIFT/image_alignment_HE_to_xen_dapi_post_xoa_v3.0.py -i ../HT185P1-S1H2L1U1.ome.tif -d ../morphology_focus_0000.ome.tif -b -s 4 -m 20