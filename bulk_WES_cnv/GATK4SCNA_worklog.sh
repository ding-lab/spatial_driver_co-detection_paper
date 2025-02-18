# bulk CNV calling steps were run one at a time.
# Running CNV calls for all HTAN samples with bulk WES involved in the Xenium variant project (tumor-only):
# Step 1 (precall) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p precall -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_allTumorBAM.list -M /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 2 (PON) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p pon -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 3 (CNV calls all using paired version) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callcn_tumor_only -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 3.5 (CNV calls normal) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callnormal -M /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 4 (plot all Tumor) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plot -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 4.5 (plot normals) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plotNormal -M /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 5 (Calls gene-level) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p geneLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 6 (Merge gene-level files to one file) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p merge -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 5.5 (Calls chr_arm-level) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 6.5 (Merge arm-level files to one file) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmmerge -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# final plot for figure 3a.
bash $austin/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh \
-p plotFinal \
-t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/HTAN_WES_geneTumorBAM.list \
-o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/HTAN/freeze_v3_2025-01-16/ \
-c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Running CNV calls for all PECGS/CPTAC samples with bulk WES involved in the Xenium variant project (tumor-only):
# Step 1 (precall) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p precall -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_allTumorBAM.list -M /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 2 (PON) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p pon -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 3 (CNV calls all using tumor-only version) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callcn -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 3.5 (CNV calls normal) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p callnormal -M /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 4 (plot all Tumor) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plot -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 4.5 (plot normals) #
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p plotNormal -M /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_normalBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/ -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 5 (Calls gene-level) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p geneLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 6 (Merge gene-level files to one file) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p merge -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 5.5 (Calls chr_arm-level) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmLevel -t /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23/PECGS_CPTAC_WES_allTumorBAM.list -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini
# Step 6.5 (Merge arm-level files to one file) ##
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/gatk_somatic.cnv.compute1.sh -p chrarmmerge -o /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/cnv/PECGS/freeze_v2_2024-07-23 -c /storage1/fs1/dinglab/Active/Projects/austins2/tools/GATK4SCNA/config/config.gatk4scna.compute1.ini

# bulk CNV calls were merged to a single file as documented in the plotting_snippets/plotting_snippets.R file.
