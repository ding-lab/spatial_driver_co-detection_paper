# running inferCNV for the snRNAseq samples that are matching the Xenium variant sample
# This worklog documents going from the the snRNA/scRNA object to the generation and attachment of inferCNV results to the metadata of the originating seurat object for plotting/analysis.
# There are some cells in the objects that are labeled as "Low quality", "unknown", or something similar. 
# For the following analysis I excluded cells with those labels from being run via inferCNV.
# This was done by assigning all neoplastic cells the label "neoplasm". 
# All normal cells were give the label "normal". All cells with the unknown label were excluded.
conda activate seurat5 
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/
bash /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/prepare_inferCNV_inputs_xen_var_v1.sh -t /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/snv_infer_cnv_inputs.txt -o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/

# The above command sets up all of the input files for each sample and structures the input directory for use in the following commands. 
# The entire input directory is then transferred from the lab server over to the WUSTL compute1 RIS HPC cluster.
# Following data transfer the following commands were run.
# This relies on the docker image: trinityctat/infercnv:1.11.2
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/run_inferCNV_compute1.trinityctat.sh -T /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-06-06/inputs/annotations_file/reference_cells.txt -D /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-06-06 -G /a.n.southard-smith/infercnv
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/run_inferCNV_compute1.trinityctat.sh -T /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-07-23/inputs/annotations_file/reference_cells_2.txt -D /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-07-23 -G /a.n.southard-smith/infercnv
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/run_inferCNV_compute1.trinityctat.sh -T /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/20241203/inputs/annotations_file/reference_cells_3.txt -D /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/20241203 -G /a.n.southard-smith/infercnv

# Wait for all of the LSF jobs to finish and check the log files to make sure that they ran successfully.
# If one of the samples fails to generate then try increasing the memory requested by the bsub call in the run_inferCNV_compute1.trinityctat.sh script: `'select[mem>240GB] rusage[mem=240GB]' -M 240GB`
# Then generate the output matrices we will rely on for plotting which cells have what copy number changes on the UMAP
# this relies on the docker image: austins2/condapython:3.9.6.base
# the output matrix table here was used to add the copy numbers that are plotted on the umap in figure 3 and the related Extended Data figure.
cd /storage1/fs1/dinglab/Active/Projects/austins2/snv_project/inferCNV/2024-06-06/outputs/matrices
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/convert-cnv-predictions-to-matrix-v3.sh -t infercnv_outs.txt -o .
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/convert-cnv-predictions-to-matrix-v3.sh -t infercnv_outs_2.txt -o .
bash /storage1/fs1/dinglab/Active/Projects/austins2/tools/inferCNV/convert-cnv-predictions-to-matrix-v3.sh -t infercnv_outs_3.txt -o .

# Transfer the complete contents of the output matrix folder from the WUSTL compute1 RIS HPC cluster back to the lab server.
conda activate seurat5
cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots
# the outputs here are where I pulled the inferCNV counts presented in the supplementary table that compares probe VAF to inferCNV copy amplification/deletion, snRNA VAF, bulk WES VAF, and bulk WES CNV call and log2ratio
bash /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/plotting_cnvs_per_cell_type.sh -t plotting_input_table.tsv -o . -g 'ATM,ARID1A,BRIP1,BRCA1,BRCA2,SMAD4,KRAS,TP53,CCND1,LCE1D,CCNE1,MCL1,ITGAE,EOMES,CD200R1,MBP,MYC,EGFR,PDGFRA,ERBB2,TCL1A,CD27,PIK3CA'
bash /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/plotting_cnvs_per_cell_type.sh -t plotting_input_table_2.tsv -o . -g 'ATM,ARID1A,BRIP1,BRCA1,BRCA2,SMAD4,KRAS,TP53,CCND1,LCE1D,CCNE1,MCL1,ITGAE,EOMES,CD200R1,MBP,MYC,EGFR,PDGFRA,ERBB2,TCL1A,CD27,PIK3CA'
bash /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/plotting_cnvs_per_cell_type.sh -t plotting_input_table_3.tsv -o . -g 'ATM,ARID1A,BRIP1,BRCA1,BRCA2,SMAD4,KRAS,TP53,CCND1,LCE1D,CCNE1,MCL1,ITGAE,EOMES,CD200R1,MBP,MYC,EGFR,PDGFRA,ERBB2,TCL1A,CD27,PIK3CA'
bash /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/plotting_cnvs_per_cell_type.sh -t plotting_input_table_4.tsv -o . -g 'PIK3CA'


