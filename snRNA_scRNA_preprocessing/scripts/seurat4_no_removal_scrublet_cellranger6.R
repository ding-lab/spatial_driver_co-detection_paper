#!/usr/bin/env Rscript --vanilla
# written by Austin Southard-Smith (a.n.southard-smith@wustl.edu)
# original script path: /diskmnt/Projects/Users/austins2/tools/seurat4_no_removal_scrublet_cellranger6.R
# scrublet version
# load required libraries
# example command: Rscript /diskmnt/Projects/Users/austins2/tools/scrublet_cleaning_seurat_v4.R -i /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/cellranger-arc2/HT224P1-S1Fc2A2N1Bma1_1/outs/raw_feature_bc_matrix -s HT224P1-S1Fc2A2N1Bma1_1 -o HT224P1-S1Fc2A2N1Bma1_1 --scrublet /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/scrublet/HT224P1-S1Fc2A2N1Bma1_1 --detailed TRUE 

library(optparse)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)

# create user options
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. cellranger output's raw matrices folder)",
              metavar="character"),
  make_option(c("--scrublet"),
              type="character",
              default=NULL,
              help="path to scrublet data folder (e.g. /diskmnt/Projects/PDX_scRNA_analysis/HTAN-BRCA/cellranger-no-intron/scrublet)",
              metavar="character"),
  make_option(c("--pre_filter"),
              type="integer",
              default=300,
              help="min number of reads per cell to prefilter",
              metavar="integer"),
  make_option(c("--nfeature_min"),
              type="integer",
              default=200,
              help="nFeature_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--nfeature_max"),
              type="integer",
              default=10000,
              help="nFeature_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--ncount_min"),
              type="integer",
              default=1000,
              help="nCount_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--ncount_max"),
              type="integer",
              default=80000,
              help="nCount_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--mito_max"),
              type="double",
              default=0.1,
              help="maximum allowed mitochondrial fraction",
              metavar="double"),
  make_option(c("-o", "--output"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-s","--sample_id"),
              type="character",
              default="single_cell_study",
              help="Name of your sample",
              metavar="character"),
  make_option(c("--pc_num"),
              type="integer",
              default=30,
              help="number of principal components to use",
              metavar="integer"),
  make_option(c("--detailed"),
              action="store_true",
              default=TRUE,
              help="relevant output plots will contain medians of nFeature_RNA, nCount_RNA, and percent.mito, and script will output a report about the number of cells processed")
              
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}
if (is.null(opt$scrublet)){
  print_help(opt_parser)
  stop("Path to scrublet is required (--scrublet).n", call.=FALSE)
}


# read in initial arguments
sample_id <- opt$sample_id
out_path <- opt$output
matrix_dir = opt$input

# make output dir if it doesn't exist
dir.create(out_path)

# get direct paths to data
barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")

# read in matrix
input <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
colnames(input) = barcode.names$V1
rownames(input) = feature.names$V2

# pre-filter and create initial seurat object
bc_sums <- colSums(input)
bg_id <- names(bc_sums[bc_sums >= opt$pre_filter])
obj = CreateSeuratObject(counts = input[,bg_id],project=opt$sample_id,min.cells = 0)


#####read in scrublet output matrix#####
scrublet.path <- paste0(opt$scrublet,"/",sample_id,"_scrublet_output_table.csv", sep="")
scrublet.df <- read.table(file = scrublet.path, sep=",", header = TRUE, row.names = "Barcodes")
obj <- AddMetaData(object = obj, metadata = scrublet.df)

### QC
# get percent mitochondrial content
obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-") / 100

# plot pre-filter metadata
#obj$percent.mito<-percent.mito
pdf(paste(out_path,"/QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
median_Features_pre_filter = median(obj$nFeature_RNA)
median_Count_pre_filter = median(obj$nCount_RNA)
median_percent.mito_pre_filter = median(obj$percent.mito)
VlnPlot(object = obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# plot metadata associations
pdf(paste0(out_path,"/FeatureScatter_in_sample_",sample_id,".pdf",sep=""),width=12,height=7)
plot1 <- FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#cell statistics
# filter step
obj<-subset(x = obj, subset = nFeature_RNA > opt$nfeature_min & nFeature_RNA < opt$nfeature_max & nCount_RNA > opt$ncount_min & nCount_RNA < opt$ncount_max & percent.mito<opt$mito_max)
obj %>% dim
num_cells_pre_doublet_removal <- length(obj$nCount_RNA)
num_cells_post_doublet_removal <- length(obj$nCount_RNA)
num_predicted_singlet_by_scrublet <- table(obj$predicted_doublet)[[1]]
num_cells_remaining_with_no_doublet_prediction <- num_cells_post_doublet_removal - num_predicted_singlet_by_scrublet


# plot post-filter metadata
pdf(paste(out_path,"/After_QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
median_Features_post_filter = median(obj$nFeature_RNA)
median_Count_post_filter = median(obj$nCount_RNA)
median_percent.mito_post_filter = median(obj$percent.mito)

VlnPlot(object = obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

dev.off()

# Run the standard workflow for visualization and clustering
obj <- SCTransform(obj, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
obj <- RunPCA(obj, npcs = opt$pc_num, verbose = FALSE)

# t-SNE and Clustering
obj <- RunUMAP(obj, reduction = "pca", dims = 1:opt$pc_num)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:opt$pc_num)
obj <- FindClusters(obj, resolution = 0.5)

# plot the clusters
pdf(paste0(out_path,"/DimPlot_",sample_id,".pdf"),useDingbats=FALSE)
DimPlot(object = obj, reduction = "umap",label=TRUE,label.size=6)
dev.off()

# save object so far
saveRDS(obj,file = paste(out_path,"/",sample_id, "_processed.rds", sep=""))
if (opt$detailed){
  sink(paste0(out_path,"/QC_statistics_",sample_id,".txt",sep=""))
  print(paste0("median nFeature_RNA pre-filtering: ", median_Features_pre_filter))
  print(paste0("median nCount_RNA pre-filtering: ", median_Count_pre_filter))
  print(paste0("median percent.mito pre-filtering: ",median_percent.mito_pre_filter))
  print(paste0("number of cells pre-doublet removal: ", num_cells_pre_doublet_removal))
  print(paste0("number of cells post-doublet removal: ", num_cells_post_doublet_removal))
  print(paste0("number of singlet cells remaining predicted by scrublet: ", num_predicted_singlet_by_scrublet))
  print(paste0("number of cells remaining without a doublet prediction (could be doublets or singlets):",num_cells_remaining_with_no_doublet_prediction))
  print(paste0("median nFeature_RNA post-filtering: ", median_Features_post_filter))
  print(paste0("median nCount_RNA post-filtering: ", median_Count_post_filter))
  print(paste0("median percent.mito post-filtering: ",median_percent.mito_post_filter))
  sink()
}


