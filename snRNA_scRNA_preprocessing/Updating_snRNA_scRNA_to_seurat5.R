#!/usr/bin/env Rscript
# conda activate seurat5
# Updating the Seurat objects that already have doublets removed.
# Used for the following samples:
# HT227P1-S1H1L1U1 
# HT242P1-S1H4L4U1 # UMAP is too spread out.
# HT284P1-S1H1A1U1 # UMAP did not generate
# SP001P1-Fp1U1
# Rscript /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/update_PDAC_objects_to_v5.0.1.R update_input_table_PDAC_v1.tsv
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table <- read.table(args[1],sep='\t',header=T,row.names="X")
# number, sample, rds_path, cell_type, cell_type_key
iterables <- rownames(input_table)
for (i in iterables) {
    sample = input_table[i,"sampleID"]
    obj_path = input_table[i,"rds_obj"]
    cell_type_input = input_table[i,"cell_type"]
    cell_type_key = input_table[i,"cell_type_key"]
    cell_type_meta <- read.table(cell_type_input, sep='\t', header=T)
    print(sample)
    obj <- readRDS(obj_path)
    DefaultAssay(obj) <- "RNA"
    obj <- DietSeurat(obj, assays = c('RNA'))
    obj <- UpdateSeuratObject(obj)
    sample_cell_type_meta <- cell_type_meta[(cell_type_meta$sample_ID == sample), ]
    tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
    print(head(tmp_df))
    obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
    obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
    obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
    all.genes <- rownames(obj)
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
    obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
    DefaultAssay(obj) <- "SCT"
    obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
    obj <- FindClusters(obj, resolution = 0.5)
    Idents(obj) <- "seurat_clusters"
    pdf(paste0("UMAP_QC_",sample,".pdf"), useDingbats = F)
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    pdf(paste0("UMAP_cell_type_xenium_varaint_v1",sample,".pdf"))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
    obj$neoplasm_normal_unknown <- NA
    obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor","PanIN","Tumor_proliferating","putative_ITPN"))] <- "neoplasm"
    obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor","PanIN","Tumor_proliferating","putative_ITPN")))] <- "normal"
    saveRDS(obj, paste0(sample,"_processed_update_no_doublet_seurat5.0.1.rds"))
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
    write.table(tmp_df, paste0(out_dir,"/",sample,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
    write.table(tmp_df, paste0(out_dir,"/",sample,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}
# HT227P1-S1H1L1U1 # Dimplot from the above was not wide enough because of the number of cell types. Replotting
# HT242P1-S1H4L4U1 # Dimplot from the above was not wide enough because of the number of cell types. Replotting
# HT284P1-S1H1A1U1 # Dimplot from the above was not wide enough because of the number of cell types. Replotting
# the above is because I forgot to do useDingbats = F
# remaking the plots below.
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1_processed_update_no_doublet_seurat5.0.1.rds")
sample = "HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
pdf(paste0("UMAP_cell_type_xenium_varaint_v1",sample,".pdf"), useDingbats = F, width=12, height=5)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT242P1-S1H1Fc2A2N1Z1_1Bmn1_processed_update_no_doublet_seurat5.0.1.rds")
sample = "HT242P1-S1H1Fc2A2N1Z1_1Bmn1"
pdf(paste0("UMAP_cell_type_xenium_varaint_v1",sample,".pdf"), useDingbats = F, width=12, height=5)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT284P1-S1H2Fc2A2N1Z1_1Bmn1_1_processed_update_no_doublet_seurat5.0.1.rds")
sample = "HT284P1-S1H2Fc2A2N1Z1_1Bmn1_1"
pdf(paste0("UMAP_cell_type_xenium_varaint_v1",sample,".pdf"), useDingbats = F, width=12, height=5)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
# doing HT270P1
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
# number, sample, rds_path, cell_type, cell_type_key
sample = "HT270P1-S1H2Fc2A2N1Bmn1_1"
obj_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/individual_no_doublets/v6/HT270P1-S1H2Fc2A2N1Bmn1_1_processed_v6.1_nod.rds"
cell_type_input = "/diskmnt/Projects/HTAN_analysis_2/PDAC/merge/v6/integration/output/no_doublets/PDAC_int_everything_no_doublets_snRNA_v6.2_metadata_20240228.tsv"
cell_type_key = "cell_type_merge_v6.2"
cell_type_meta <- read.table(cell_type_input, sep='\t', header=T)
print(sample)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
obj <- UpdateSeuratObject(obj)
sample_cell_type_meta <- cell_type_meta[(cell_type_meta$sample_ID == sample), ]
tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
print(head(tmp_df))
obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
all.genes <- rownames(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
pdf(paste0("UMAP_QC_",sample,".pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
pdf(paste0("UMAP_cell_type_xenium_varaint_v1",sample,".pdf"), useDingbats = F, width=12, height=5)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor","PanIN","Tumor_proliferating","putative_ITPN"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor","PanIN","Tumor_proliferating","putative_ITPN")))] <- "normal"
saveRDS(obj, paste0(sample,"_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
# generating the single cell seurat object for the mCRC case that needs to be subset from the pancan-snATAC matchign snRNA object
# conda activate seurat5
#!/usr/bin/env Rscript
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
sample = "HT260C1-Th1K1Fc2A2N1Bma1_1"
piece = "HT260C1-Th1K1"
obj_path = "/diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/CRC/merge-final/v7/CRC_merge_obj_no_doublets_v7.rds"
cell_type_input = "/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snRNA/All_snRNA_samples_metadata_data_freeze_v7.0.tsv"
cell_type_key = "cell_type.harmonized.cancer"
cell_type_meta <- read.table(cell_type_input, sep='\t', header=T)
print(sample)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
obj <- UpdateSeuratObject(obj)
obj <- subset(obj, subset = Piece_ID == piece)
sample_cell_type_meta <- cell_type_meta[cell_type_meta$Sample_RNA == sample, ]
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
head(tmp_df)
obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
obj <- RenameCells(obj, new.names = paste0(sample,"_",obj$gex_barcode))
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
all.genes <- rownames(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
pdf(paste0(out_dir,"UMAP_QC_",sample,".pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_rna", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_rna", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_atac", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_atac", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample,".pdf"))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
#                  0   1   2   3   4   5   6   7   8   9  10
# B-cells          0   0   0   0   0   1   0   0   0   0   0
# Cholangiocytes   0   0   0   0   0   0   0   0   8   0   0
# Endothelial      0   0   0   0   0   2   0   0   4   0  14
# Fibroblasts      0   0   0   0   0   6   0   0  38   0   2
# Hepatocytes      0   0   0   0   0   1   0   0   3   0   0
# Low quality      0   0   0   1   0   1   1   1   2   0   1
# Macrophages      2   0   0   0   0  32   0   0   6   1   2            
# Plasma           0   0   1   0   0   2   0   1   1   0   1
# T-cells          1   0   0   0   0  17   0   0   5   0   0
# Tumor          628 620 506 437 366 143 203 145  70 101  23
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal"
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Low quality"))] <- "unknown"
saveRDS(obj, paste0(out_dir,sample,"_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
# this object has been stringently filtered and as a result there are very few tumor cells present in it. Clara has provided an alternative.
#!/usr/bin/env Rscript
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
sample = "HT260C1-Th1K1Fc2A2N1Bma1_1"
piece = "HT260C1-Th1K1"
obj_path = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/28-snRNA/28_0-snRNA_obj/HT260C1-Th1K1/HT260C1-Th1K1_processed.rds"
cell_type_input = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/28-snRNA/28_4-sample_level_cluster/metadata_table_HT260C1-Th1K1.tsv"
cell_type_key = "cell_type_final"
cell_type_meta <- read.table(cell_type_input, sep='\t', header=T)
print(sample)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
obj <- UpdateSeuratObject(obj)
scrublet.df <- read.table("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/28-snRNA/28_2-snRNA_scrublet/output/cellranger/HT260C1-Th1K1/HT260C1-Th1K1_scrublet_output_table.csv", sep=",", header = TRUE, row.names = "Barcodes")
obj <- AddMetaData(object = obj, metadata = scrublet.df)
doublet_finder_table <- read.table("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/28-snRNA/28_3-snRNA_DoubletFinder/HT260C1-Th1K1/HT260C1-Th1K1_DoubletFinder.tsv", sep="\t", header = TRUE, row.names = "Barcodes")
obj <- AddMetaData(object = obj, metadata = doublet_finder_table)
obj <- subset(x = obj, subset = predicted_doublet == 'False'| is.na(obj$predicted_doublet))
obj <- subset(obj, subset = DF.classifications_0.25_0.09_803 == "Singlet")
sample_cell_type_meta <- cell_type_meta
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
head(tmp_df)
obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
obj <- RenameCells(obj, new.names = paste0(sample,"_",obj$original_RNA_barcode))
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
all.genes <- rownames(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
pdf(paste0(out_dir,"/","UMAP_QC_",sample,"_clara.pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "pANN_0.25_0.09_803", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "DF.classifications_0.25_0.09_803", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","UMAP_cell_type_xenium_varaint_v1",sample,"_clara.pdf"))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
#                  0    1    2    3    4    5    6    7    8    9   10   11
# B                2    1    1    2    0    1    0    0    0    0    0  102
# Cholangiocyte    0    0    0    0    0    0    0    0    2    0    0    0
# Endothelial      0    2    0    0    0    0    0    0    1    0  225    0
# Fibroblast       0    0    0    0    0    0    0    0  394    0    0    0
# Hepatocyte       0    0    0    0    0    1    0    0    0    0    0    1
# Macrophage       1    1    3  920    0    3    0    0    0    0    0    8
# T                2    3    2    0    3  806    0    0    0    0    0    5
# Tumor         1363 1342 1029    2  815    2  697  414    1  257    1   91
#                 12   13
# B                0    0
# Cholangiocyte    0   17
# Endothelial      0    0
# Fibroblast       2    0
# Hepatocyte       1   59
# Macrophage       5    0
# T                8    0
# Tumor          129    2
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal"
#obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Low quality"))] <- "unknown"
saveRDS(obj, paste0(out_dir,"/",sample,"_clara_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample,"_clara_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample,"_clara_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
#The UMAP of this one looks better and there is a much higher diversity of cell types here. Will use the object provided by Clara instead.

# generating the single cell objects for the HTAN BRCA cases:
# the breast cancer data is split out across 3 objects. 
# 1 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_multiome_snRNA_v10312023_merged_obj_varFeatOnly_annotated.rds
# 2 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_snRNA_merged_obj_celltypeannotation.rds
# 3 /diskmnt/Projects/HTAN_analysis_2/PDAC/HTAN_BRCA_scRNA_merged_obj_v2.rds
# 
# 1 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_multiome_snRNA_v10312023_merged_obj_varFeatOnly_annotated_metadata.tsv
# 2 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_snRNA_merged_obj_celltypeannotation_metadata.tsv
# 3 /diskmnt/Projects/HTAN_analysis_2/PDAC/HTAN_BRCA_scRNA_merged_obj__transfered_old_metadata_to_new_object_metadata_20230913.tsv
#!/usr/bin/env Rscript
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
sample = ""
piece = ""
cell_type_key = "cell_type_rj"
obj_1_path = "/diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_multiome_snRNA_v10312023_merged_obj_varFeatOnly_annotated.rds"
meta_1_path = "/diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_multiome_snRNA_v10312023_merged_obj_varFeatOnly_annotated_metadata.tsv"
meta_1 <- read.table(meta_1_path, sep='\t', header=T)
unique(meta_1$Sample_ID)
# [1] NA                                 "HT110B1-1693N1K2Y3_1N1Z1_1Bmn1_1"
# [3] "HT235B1-S1H1Fc2A2N1Z1_1Bmn1_1"    "HT243B1-S1H4Fc2A2N1Z1_1Bmn1"
# [5] "HT263B1-S1H1A3N1Z2_1Bmn1_1"       "HT265B1-S1H1Fc2A2N1Bmn1_1"
# [7] "HT268B1-Th1H3Fc2A2N1Z1_1Bmn1_1"   "HT271B1-S1H3Fc2A5N1Bmn1_1"
# [9] "HT297B1-S1H1Fc2A2N1Z1Bmn1_1"      "HT305B1-S1H1Fc2A2_1N1Bmn1_1"
# [11] "HT308B1-S1V1Fc2A2N1Z1_1Bmn1_1"    "HT309B1-S1H3A3Y1Nd1_1Z1_1Bmn1_1"
# [13] "HT323B1-S1H1Fc2A2N1Z1_1Bmn1_1"    "HT339B1-S1H3Fc2A2N1Z1_1Bmn1_1"
# [15] "HT339B2-S1H2A3Y1Nd1_1Z1_1Bmn1_1"  "HT355B1-S1H1A3N1Z1_1Bmn1_1"
# [17] "HT365B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT365B1-S1H1Fc2A2N1Z1_1Bmn1_1"
# [19] "HT372B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT374B1-S1H3Fc2A2_1N1Z1_1Bmn1_1"
# [21] "HT378B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT378B1-S1H2A3Y1N1Z1_1Bmn1_1"
# [23] "HT384B1-M1_1N1Z1_1Bmn1_1"         "HT384B1-S1H1A3_1N1Z1_1Bmn1_1"
# [25] "HT384B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT389B1-M1A3N1Z1_1Bmn1_1"
# [27] "HT397B1-S1H4A4Y1N1Z1_1Bmn1_1"     "HT423B1-S1H2A2Y1N1Z1_1Bmn1_1"
# [29] "HT425B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT461B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
# [31] "HT480B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT482B1-S1H1A3Y1N1Z1_1Bmn1_1"
# [33] "HT486B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT497B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
# [35] "HT514B1-S1H3A3Y1Nd1_1Z1_1Bmn1_1"  "HT517B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
# [37] "HT519B1-S1H3A3Y1Nd1_1Z1_1Bmn1_1"  "HT533B1-H1Fc1Nd1_1Z1_1Bmn1_1"
# [39] "HT545B1-H1Fc2Nd1_1Z1_1Bmn1_1"     "HT565B1-S1H6A3Y1Nd1_1Z1_1Bmn1_1"
# [41] "HT577B1-H1Fc2Nd1_1Z1_1Bmn1_1"     "HT591B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
# [43] "HT617B1-S1H4A3Y1Nd1_1Z1_1Bmn1_1"  "HT632B1-S1H4A3Y1Nd1_1Z1_1Bmn1_1"
# [45] "HT648B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"  "HT686B1-S1H4A3Y1Nd1_1Z1_1Bmn1_1"
# [47] "HT701B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
unique(meta_1$orig.ident)
# [1] "BC003B2-S1Fc1Y1Nd1_1Z1_1Bmn1_1"   "HT110B1-1693N1K2Y3_1N1Z1_1Bmn1_1"
# [3] "HT235B1-S1H1Fc2A2N1Z1_1Bmn1_1"    "HT243B1-S1H4Fc2A2N1Z1_1Bmn1"     
# [5] "HT263B1-S1H1A3N1Z2_1Bmn1_1"       "HT265B1-S1H1Fc2A2N1Bmn1_1"       
# [7] "HT268B1-Th1H3Fc2A2N1Z1_1Bmn1_1"   "HT271B1-S1H3Fc2A5N1Bmn1_1"       
# [9] "HT297B1-S1H1Fc2A2N1Z1Bmn1_1"      "HT305B1-S1H1Fc2A2_1N1Bmn1_1"     
# [11] "HT308B1-S1V1Fc2A2N1Z1_1Bmn1_1"    "HT309B1-S1H3A3Y1Nd1_1Z1_1Bmn1_1" 
# [13] "HT323B1-S1H1Fc2A2N1Z1_1Bmn1_1"    "HT339B1-S1H3Fc2A2N1Z1_1Bmn1_1"   
# [15] "HT339B2-S1H2A3Y1Nd1_1Z1_1Bmn1_1"  "HT355B1-S1H1A3N1Z1_1Bmn1_1"      
# [17] "HT365B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT365B1-S1H1Fc2A2N1Z1_1Bmn1_1"   
# [19] "HT372B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT374B1-S1H3Fc2A2_1N1Z1_1Bmn1_1" 
# [21] "HT378B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT378B1-S1H2A3Y1N1Z1_1Bmn1_1"    
# [23] "HT384B1-M1_1N1Z1_1Bmn1_1"         "HT384B1-S1H1A3_1N1Z1_1Bmn1_1"    
# [25] "HT384B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT389B1-M1A3N1Z1_1Bmn1_1"        
# [27] "HT397B1-S1H4A4Y1N1Z1_1Bmn1_1"     "HT423B1-S1H2A2Y1N1Z1_1Bmn1_1"    
# [29] "HT425B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT461B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1" 
# [31] "HT480B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT482B1-S1H1A3Y1N1Z1_1Bmn1_1"    
# [33] "HT486B1-S1H1A3Y1N1Z1_1Bmn1_1"     "HT497B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1" 
# [35] "HT514B1-S1H3A3Y1Nd1_1Z1_1Bmn1_1"  "HT517B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1" 
# [37] "HT519B1-S1H3A3Y1Nd1_1Z1_1Bmn1_1"  "HT533B1-H1Fc1Nd1_1Z1_1Bmn1_1"    
# [39] "HT545B1-H1Fc2Nd1_1Z1_1Bmn1_1"     "HT565B1-S1H6A3Y1Nd1_1Z1_1Bmn1_1" 
# [41] "HT577B1-H1Fc2Nd1_1Z1_1Bmn1_1"     "HT591B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1" 
# [43] "HT617B1-S1H4A3Y1Nd1_1Z1_1Bmn1_1"  "HT632B1-S1H4A3Y1Nd1_1Z1_1Bmn1_1" 
# [45] "HT648B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"  "HT686B1-S1H4A3Y1Nd1_1Z1_1Bmn1_1" 
# [47] "HT701B1-S1H1A3Y1Nd1_1Z1_1Bmn1_1"
obj_2_path = "/diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_snRNA_merged_obj_celltypeannotation.rds"
meta_2_path = "/diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_snRNA_merged_obj_celltypeannotation_metadata.tsv"
meta_2 <- read.table(meta_2_path, sep='\t', header=T)
unique(meta_2$Sample_ID)
# [1] NA                                "HT-023-B-Slice1_3v3sn-lib1"     
# [3] "HT027B1-S1R1A4G1"                "HT-027-B-PunchA1_3v3sn-lib1"    
# [5] "HT-027-B-PunchB1_3v3sn-lib1"     "HT-027-B-PunchD1_3v3sn-lib1"    
# [7] "HT-027-B-Slice1fresh_snRNA-lib1" "HT-027-B-Slice1_snRNA-lib1"     
# [9] "HT029B1-PCH1Y1_XBn1_1"           "HT035B-XBn1_1"                  
# [11] "HT036B1-S1PEH2YN1Z1B1"           "HT036B1-S2PGH1Y1_XBn1_1"        
# [13] "HT065B1-S1H1A3N1"                "HT088B1-S1H1A2K2Y3_XBn1_1"      
# [15] "HT088B1-S1H2A2Y1"                "HT105B1-1679N1K2Y1_1N1Z1_1B1"   
# [17] "HT110B1-S1H4A3Y2_XBn1_1"         "HT128B1-S1H3A2K2Y1-XBn1"        
# [19] "HT128B1-S1H4A2K2Y1_1N1Z1_1B1"    "HT137B1-S1H7A2K2Y1_XBn1_1"      
# [21] "HT141B1-S1H1A3Y2N1Z1_1Bn1"       "HT163B1-S1H6A3-XBn1"            
# [23] "HT185B1-S1H2A2K1G1Z1_1B1"        "HT185B1-S1H3A2K1G1Z1_1B1"       
# [25] "HT185B1-S1H6A2K1G1Z1_1B1"        "HT206B1-S1H4A2Y1_1N1Z1_1B1"     
# [27] "HT206B1-S1H4A2Y1_1N1Z1_1Bn1_1"   "HT206B1-S1H4A2Y2_1N3Z1_1Bn1"    
# [29] "HT214B1-S1H2A2Y1_1N1Z1_1B1"      "HT214B1-S1H2A3Y1_3N1Z1_1Bn2_1"  
# [31] "HT217B1-S1H1A2Y1_1N1Z1_1B1"      "HT217B1-S1H1A3Y1_1N1Z1_1Bn1_1"  
# [33] "HT235B1-H1S1Fc1_1N1Z1Bn1_1"      "HT243B1-H3A2Fc1N1Z1_1Bn1_1"     
# [35] "HT262B1-S1H3A3N1Z1_2Bn2"         "HTAN_1408-06-3x3sn-lib1"
unique(meta_2$orig.ident)
# [1] "HT-016-B-S1_R2_2_snRNA-lib1"     "HT-023-B-Slice1_3v3sn-lib1"     
# [3] "HT027B1-S1R1A4G1"                "HT-027-B-PunchA1_3v3sn-lib1"    
# [5] "HT-027-B-PunchB1_3v3sn-lib1"     "HT-027-B-PunchD1_3v3sn-lib1"    
# [7] "HT-027-B-Slice1fresh_snRNA-lib1" "HT-027-B-Slice1_snRNA-lib1"     
# [9] "HT029B1-PCH1Y1_XBn1_1"           "HT035B-XBn1_1"                  
# [11] "HT036B1-S1PEH2YN1Z1B1"           "HT036B1-S2PGH1Y1_XBn1_1"        
# [13] "HT065B1-S1H1A3N1"                "HT088B1-S1H1A2K2Y3_XBn1_1"      
# [15] "HT088B1-S1H2A2Y1"                "HT105B1-1679N1K2Y1_1N1Z1_1B1"   
# [17] "HT110B1-S1H4A3Y2_XBn1_1"         "HT128B1-S1H3A2K2Y1-XBn1"        
# [19] "HT128B1-S1H4A2K2Y1_1N1Z1_1B1"    "HT137B1-S1H7A2K2Y1_XBn1_1"      
# [21] "HT141B1-S1H1A3Y2N1Z1_1Bn1"       "HT163B1-S1H6A3-XBn1"            
# [23] "HT185B1-S1H2A2K1G1Z1_1B1"        "HT185B1-S1H3A2K1G1Z1_1B1"       
# [25] "HT185B1-S1H6A2K1G1Z1_1B1"        "HT206B1-S1H4A2Y1_1N1Z1_1B1"     
# [27] "HT206B1-S1H4A2Y1_1N1Z1_1Bn1_1"   "HT206B1-S1H4A2Y2_1N3Z1_1Bn1"    
# [29] "HT214B1-S1H2A2Y1_1N1Z1_1B1"      "HT214B1-S1H2A3Y1_3N1Z1_1Bn2_1"  
# [31] "HT217B1-S1H1A2Y1_1N1Z1_1B1"      "HT217B1-S1H1A3Y1_1N1Z1_1Bn1_1"  
# [33] "HT235B1-H1S1Fc1_1N1Z1Bn1_1"      "HT243B1-H3A2Fc1N1Z1_1Bn1_1"     
# [35] "HT262B1-S1H3A3N1Z1_2Bn2"         "HTAN_1408-06-3x3sn-lib1"
obj_3_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/HTAN_BRCA_scRNA_merged_obj_v2.rds"
meta_3_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/HTAN_BRCA_scRNA_merged_obj__transfered_old_metadata_to_new_object_metadata_20230913.tsv"
meta_3 <- read.table(meta_3_path, sep='\t', header=T)
unique(meta_3$Sample_ID)
# [1] NA                             "HT062B1-S1PBA1A1Z1B1"        
# [3] "HT062B1-S1R1A1Z1B1"           "HT065B1-S1H1A2A1Z1B1"        
# [5] "HT065B1-S1H7A2A1Z1B1"         "HT067B1-S1H2A2A1Z1B1"        
# [7] "HT067B1-S1H5A2A1Z1B1"         "HT068B1-S1H1A2A1Z1B1"        
# [9] "HT068B1-S1H2A2A1Z1B1"         "HT069B1-S1H3A2A1G1Z1B1"      
# [11] "HT069B1-S1H9A2A1G1Z1B1"       "HT077B1-S1H1A2K1G1Z1B1"      
# [13] "HT077B1-S1H3A2K1G1Z1B1"       "HT077B1-S1H7A2K1G1Z1B1"      
# [15] "HT084B1-S1H8A2A1G1Z1B1"       "HT088B1-S1H1A2K1G1Z1_1B2"    
# [17] "HT088B1-S1H4A2K1G1Z1_1B2"     "HT103B1-S1H4A2K1G1Z1B1"      
# [19] "HT103B1-S1H7A2K1G1Z1B1"       "HT105B1-1679N1K1G1Z1B1"      
# [21] "HT105B1-S1H1A2K1G1Z1B1"       "HT105B1-S1H3A2K1G1Z1B1"      
# [23] "HT110B1-N1K1Z1B1"             "HT110B1-S1H1A2K1G1Z1B1"      
# [25] "HT110B1-S1H4A2K1G1Z1B1"       "HT144B1-S1H9A2K1Z1B1"        
# [27] "HT154B1-S1H1A2K1G1Z1B1"       "HT154B1-S1H3A2K1G1Z1B1"      
# [29] "HT154B1-S1H5A2K1G1Z1B1"       "HT163B1-S1H2A2K1G1Z1B1"      
# [31] "HT163B1-S1H6A2K1G1Z1B1"       "HT171B1-S1H1A2K1G1Z1_1B1"    
# [33] "HT171B1-S1H8A2K1G1Z1_1B1"     "HT265B1-S1H4A1K2G1Z1_1Bc1_1" 
# [35] "HT268B1-TH1M1A2K2G1Z1_1Bc1_1" "HT268B1-TH2H1A4K2G1Z1_1Bc1_1"
unique(meta_3$orig.ident)
# [1] "HT062B1-S1PAA1A1Z1B1"         "HT062B1-S1PBA1A1Z1B1"        
# [3] "HT062B1-S1R1A1Z1B1"           "HT065B1-S1H1A2A1Z1B1"        
# [5] "HT065B1-S1H7A2A1Z1B1"         "HT067B1-S1H2A2A1Z1B1"        
# [7] "HT067B1-S1H5A2A1Z1B1"         "HT068B1-S1H1A2A1Z1B1"        
# [9] "HT068B1-S1H2A2A1Z1B1"         "HT069B1-S1H3A2A1G1Z1B1"      
# [11] "HT069B1-S1H9A2A1G1Z1B1"       "HT077B1-S1H1A2K1G1Z1B1"      
# [13] "HT077B1-S1H3A2K1G1Z1B1"       "HT077B1-S1H7A2K1G1Z1B1"      
# [15] "HT084B1-S1H8A2A1G1Z1B1"       "HT088B1-S1H1A2K1G1Z1_1B2"    
# [17] "HT088B1-S1H4A2K1G1Z1_1B2"     "HT103B1-S1H4A2K1G1Z1B1"      
# [19] "HT103B1-S1H7A2K1G1Z1B1"       "HT105B1-1679N1K1G1Z1B1"      
# [21] "HT105B1-S1H1A2K1G1Z1B1"       "HT105B1-S1H3A2K1G1Z1B1"      
# [23] "HT110B1-N1K1Z1B1"             "HT110B1-S1H1A2K1G1Z1B1"      
# [25] "HT110B1-S1H4A2K1G1Z1B1"       "HT144B1-S1H9A2K1Z1B1"        
# [27] "HT154B1-S1H1A2K1G1Z1B1"       "HT154B1-S1H3A2K1G1Z1B1"      
# [29] "HT154B1-S1H5A2K1G1Z1B1"       "HT163B1-S1H2A2K1G1Z1B1"      
# [31] "HT163B1-S1H6A2K1G1Z1B1"       "HT171B1-S1H1A2K1G1Z1_1B1"    
# [33] "HT171B1-S1H8A2K1G1Z1_1B1"     "HT265B1-S1H4A1K2G1Z1_1Bc1_1" 
# [35] "HT268B1-TH1M1A2K2G1Z1_1Bc1_1" "HT268B1-TH2H1A4K2G1Z1_1Bc1_1"
# there are some samples that have missing sample_ID fields so you need to rely on the orig.ident label added on at the individual object stage.
samples_to_look_up <- c('HT243B','HT268B','HT305B','HT308B','HT425B','HT065B')
meta_1_samples <- c()
for (sample_id in samples_to_look_up) {
    meta_1_samples <- c(meta_1_samples, unique(meta_1$orig.ident)[grepl(sample_id,unique(meta_1$orig.ident))])
}
print(meta_1_samples)
# [1] "HT243B1-S1H4Fc2A2N1Z1_1Bmn1"    "HT268B1-Th1H3Fc2A2N1Z1_1Bmn1_1"
# [3] "HT305B1-S1H1Fc2A2_1N1Bmn1_1"    "HT308B1-S1V1Fc2A2N1Z1_1Bmn1_1"
# [5] "HT425B1-S1H1A3Y1N1Z1_1Bmn1_1"
meta_2_samples <- c()
for (sample_id in samples_to_look_up) {
    meta_2_samples <- c(meta_2_samples, unique(meta_2$orig.ident)[grepl(sample_id,unique(meta_2$orig.ident))])
}
print(meta_2_samples)
# [1] "HT243B1-H3A2Fc1N1Z1_1Bn1_1" "HT065B1-S1H1A3N1" 
meta_3_samples <- c()
for (sample_id in samples_to_look_up) {
    meta_3_samples <- c(meta_3_samples, unique(meta_3$orig.ident)[grepl(sample_id,unique(meta_3$orig.ident))])
}
print(meta_3_samples)
# [1] "HT268B1-TH1M1A2K2G1Z1_1Bc1_1" "HT268B1-TH2H1A4K2G1Z1_1Bc1_1"
# [3] "HT065B1-S1H1A2A1Z1B1"         "HT065B1-S1H7A2A1Z1B1"
#These are the Xenium sample_IDs
# HT243B1-S1H1A4U1
# HT268B1-Th1H3L1U1
# HT305B1-S1H5A1U1
# HT308B1-S1H5A4U1
# HT425B1-S1H1Fp1U1
# HT065B1-H1A1A4U1
meta_1_sample_to_include <- c('HT243B1-S1H4Fc2A2N1Z1_1Bmn1','HT268B1-Th1H3Fc2A2N1Z1_1Bmn1_1',"HT305B1-S1H1Fc2A2_1N1Bmn1_1","HT308B1-S1V1Fc2A2N1Z1_1Bmn1_1","HT425B1-S1H1A3Y1N1Z1_1Bmn1_1")
meta_2_sample_to_include <- c('HT243B1-H3A2Fc1N1Z1_1Bn1_1',"HT065B1-S1H1A3N1")
meta_3_sample_to_include <- c("HT065B1-S1H1A2A1Z1B1","HT065B1-S1H7A2A1Z1B1")
obj_1 <- readRDS(obj_1_path)
obj_1 <- UpdateSeuratObject(obj_1)
obj <- DietSeurat(obj, assays = c('RNA'))
cell_type_key = "cell_type_rj"
DefaultAssay(obj_1) <- "RNA"
for (sample.id in meta_1_sample_to_include) {
    print(paste0("subsetting ", sample.id))
    obj <- subset(obj_1, subset = orig.ident == sample.id)
    obj <- subset(x = obj, subset = predicted_doublet == 'False'| is.na(obj$predicted_doublet))
    obj <- DietSeurat(obj, assays = c('RNA'))
    sample_cell_type_meta <- meta_1[meta_1$orig.ident == sample.id, ]
    tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
    head(tmp_df)
    obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
    obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
    obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
    obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
    obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
    obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
    all.genes <- rownames(obj)
    print("starting normalization")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
    obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
    DefaultAssay(obj) <- "SCT"
    obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
    obj <- FindClusters(obj, resolution = 0.5)
    Idents(obj) <- "seurat_clusters"
    print("plotting world domination")
    pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
    obj$neoplasm_normal_unknown <- NA
    obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
    obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal"
    saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
    write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
    write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}
obj_2 <- readRDS(obj_2_path)
obj_2 <- UpdateSeuratObject(obj_2)
cell_type_key = "cell_type_rj"
DefaultAssay(obj_2) <- "RNA"
for (sample.id in meta_2_sample_to_include) {
    print(paste0("subsetting ", sample.id))
    obj <- subset(obj_2, subset = orig.ident == sample.id)
    obj <- subset(x = obj, subset = predicted_doublet == 'False'| is.na(obj$predicted_doublet))
    obj <- DietSeurat(obj, assays = c('RNA'))
    sample_cell_type_meta <- meta_2[meta_2$orig.ident == sample.id, ]
    tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
    head(tmp_df)
    obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
    obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
    obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
    obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
    obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
    obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
    all.genes <- rownames(obj)
    print("starting normalization")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
    obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
    DefaultAssay(obj) <- "SCT"
    obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
    obj <- FindClusters(obj, resolution = 0.5)
    Idents(obj) <- "seurat_clusters"
    print("plotting world domination")
    pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
    obj$neoplasm_normal_unknown <- NA
    obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
    obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal"
    saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
    write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
    write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}
obj_3 <- readRDS(obj_3_path)
obj_3 <- UpdateSeuratObject(obj_3)
cell_type_key = "cell_type_specific"
DefaultAssay(obj_3) <- "RNA"
for (sample.id in meta_3_sample_to_include) {
    print(paste0("subsetting ", sample.id))
    obj <- subset(obj_3, subset = orig.ident == sample.id)
    obj <- subset(x = obj, subset = predicted_doublet == 'False'| is.na(obj$predicted_doublet))
    obj <- DietSeurat(obj, assays = c('RNA'))
    sample_cell_type_meta <- meta_3[meta_3$orig.ident == sample.id, ]
    tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = rownames(sample_cell_type_meta))
    head(tmp_df)
    obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
    obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
    obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
    obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
    obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
    obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
    all.genes <- rownames(obj)
    print("starting normalization")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
    obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
    DefaultAssay(obj) <- "SCT"
    obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
    obj <- FindClusters(obj, resolution = 0.5)
    Idents(obj) <- "seurat_clusters"
    print("plotting world domination")
    pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
    print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
    obj$neoplasm_normal_unknown <- NA
    obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
    obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal"
    saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
    write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
    write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}
path_list <- c(
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A2A1Z1B1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H7A2A1Z1B1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT242P1-S1H1Fc2A2N1Z1_1Bmn1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT243B1-H3A2Fc1N1Z1_1Bn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT243B1-S1H4Fc2A2N1Z1_1Bmn1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT260C1-Th1K1Fc2A2N1Bma1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT260C1-Th1K1Fc2A2N1Bma1_1_clara_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT268B1-Th1H3Fc2A2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT284P1-S1H2Fc2A2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT305B1-S1H1Fc2A2_1N1Bmn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT308B1-S1V1Fc2A2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv",
    "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/PM1467P1-T1Y2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
)
for (path in path_list) {
    cell_types <- read.table(path, sep='\t',header=F)
    print(path)
    print(unique(cell_types$V2))
}
# the cell types in each sample:
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A2A1Z1B1_cell_type_xenium_variant_v1.tsv"
# [1] "B"                  "Plasma"             "Tumor"             
# [4] "Macrophage"         "CD4_T"              "mCAF"              
# [7] "vCAF"               "Endothelial"        "Luminal_mature"    
# [10] "pDC"                "Monocyte"           "NK"                
# [13] "CD4_CD8_T"          "Basal_progenitor"   NA                  
# [16] "Luminal_progenitor" "Treg"               "dCAF"              
# [19] "Mast"               "cDC2"               "Proliferating_T"   
# [22] "CD8_T_cytotoxic"    "cDC1"               "cCAF"              
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"               "T_cell"              "Fibroblast"         
# [4] "Mono_Macro"          "Plasma"              "Adipocyte_CAF_mixed"
# [7] "Basal_myoepithelial" "B_cell"              "Luminal_progenitor" 
# [10] "Endothelial_cells"   "Mast?"               "NK"                 
# [13] "Mono_Macro?"        
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H7A2A1Z1B1_cell_type_xenium_variant_v1.tsv"
# [1] "vCAF"               "mCAF"               "Endothelial"       
# [4] "Tumor"              "Luminal_progenitor" "Macrophage"        
# [7] "Luminal_mature"     "Basal_progenitor"   "CD4_CD8_T"         
# [10] "Monocyte"           "Treg"               "cDC2"              
# [13] "Plasma"             "dCAF"               "NK"                
# [16] "B"                  "CD8_T_cytotoxic"    NA                  
# [19] "CD4_T"              "Proliferating_T"    "pDC"               
# [22] "Mast"               "cDC1"               "cCAF"              
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                  "Islet_alpha"            "Vascular_smooth_muscle"
# [4] "B_cell"                 "Duct_like_1"            "iCAF"                  
# [7] "cDC2"                   "apCAF"                  "Pericyte"              
# [10] "TAM"                    "Tuft"                   "Tumor_proliferating"   
# [13] "myCAF"                  "Islet_beta"             "Treg"                  
# [16] "CXCR4+CAF"              "Duct_like_2"            "PanIN"                 
# [19] "pDC"                    "Glial"                  "Proliferating_myeloid" 
# [22] "CD4+T"                  "AXL_DC"                 "CD8+T"                 
# [25] "putative_ADM"           "Plasma"                 "NK"                    
# [28] "Islet_delta_epsilon"    "preAP_mesothelial_CAF"  "Endothelial"           
# [31] "mregDC"                 "cDC1"                   "Islet_gamma"           
# [34] "fat_metabolism_CAF"     "Monocyte"               "ILC"                   
# [37] "Mast"                   "Lymphatic"              "Ampullary_duct_2"      
# [40] "Proliferating_T_cells"  "Acinar"                 "Ampullary_duct_1"      
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT242P1-S1H1Fc2A2N1Z1_1Bmn1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                  "Plasma"                 "myCAF"                 
# [4] "CD8+T"                  "NK"                     "Proliferating_T_cells" 
# [7] "TAM"                    "Treg"                   "PanIN"                 
# [10] "Islet_alpha"            "Islet_beta"             "CD4+T"                 
# [13] "CXCR4+CAF"              "Tumor_proliferating"    "Mast"                  
# [16] "preAP_mesothelial_CAF"  "Islet_delta_epsilon"    "cDC2"                  
# [19] "Monocyte"               "AXL_DC"                 "apCAF"                 
# [22] "iCAF"                   "putative_ADM"           "Basophil"              
# [25] "Duct_like_1"            "mregDC"                 "B_cell"                
# [28] "Duct_like_2"            "Islet_gamma"            "cDC1"                  
# [31] "Proliferating_myeloid"  "Endothelial"            "Tuft"                  
# [34] "Acinar"                 "Glial"                  "Pericyte"              
# [37] "pDC"                    "Islet_beta_delta"       "Lymphatic"             
# [40] "Ampullary_duct_1"       "Vascular_smooth_muscle" "ILC"                   
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT243B1-H3A2Fc1N1Z1_1Bn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"               "Basal_myoepithelial" "Mono_Macro"         
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT243B1-S1H4Fc2A2N1Z1_1Bmn1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                        "CAF"                         
# [3] "Myeloid"                      "Epithelial (Lum Progenitor?)"
# [5] "T_NK"                         "B"                           
# [7] "Basal_Myoepithelial"          "Plasma"                      
# [9] "CAF_Perivascular"             "Endothelial"                 
# [11] "Mast"                         "CAF_adipocyte?"              
# [13] "pDC"                          "Epithelial (Tumor?)"         
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT260C1-Th1K1Fc2A2N1Bma1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"          "Hepatocytes"    "Fibroblasts"    "Macrophages"   
# [5] "Endothelial"    "T-cells"        "Plasma"         "Low quality"   
# [9] "Cholangiocytes" "B-cells"       
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT260C1-Th1K1Fc2A2N1Bma1_1_clara_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"         "T"             "Macrophage"    "Fibroblast"   
# [5] "B"             "Endothelial"   "Cholangiocyte" "Hepatocyte"   
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT268B1-Th1H3Fc2A2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"            "T_NK"             "Myeloid"          "pDC"             
# [5] "Endothelial"      "CAF"              "CAF_Perivascular" "Plasma"          
# [9] "Mast"             "B"               
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT284P1-S1H2Fc2A2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Vascular_smooth_muscle" "Islet_alpha"            "fat_metabolism_CAF"    
# [4] "myCAF"                  "Treg"                   "Islet_beta"            
# [7] "CXCR4+CAF"              "TAM"                    "Islet_beta_delta"      
# [10] "putative_ITPN"          "Lymphatic"              "Duct_like_1"           
# [13] "Duct_like_2"            "Tumor"                  "iCAF"                  
# [16] "cDC2"                   "Glial"                  "CD4+T"                 
# [19] "Endothelial"            "Ampullary_duct_1"       "Islet_delta_epsilon"   
# [22] "CD8+T"                  "NK"                     "cDC1"                  
# [25] "ILC"                    "AXL_DC"                 "Pericyte"              
# [28] "Acinar"                 "putative_ADM"           "B_cell"                
# [31] "PanIN"                  "Ampullary_duct_2"       "Mast"                  
# [34] "Plasma"                 "apCAF"                  "Islet_gamma"           
# [37] "mregDC"                 "Adipocyte"              "Tumor_proliferating"   
# [40] "Proliferating_myeloid"  "pDC"                    "Basophil"              
# [43] "preAP_mesothelial_CAF" 
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT305B1-S1H1Fc2A2_1N1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                        "CAF"                         
# [3] "Epithelial (Lum Progenitor?)" "Basal_Myoepithelial"         
# [5] "Endothelial"                  "Myeloid"                     
# [7] "CAF_Perivascular"             "T_NK"                        
# [9] "Unknown_1"                    "Epithelial (Tumor?)"         
# [11] "Mast"                         "CAF_adipocyte?"              
# [13] "B"                            "pDC"                         
# [15] "Plasma"                      
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT308B1-S1V1Fc2A2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                        "CAF"                         
# [3] "Myeloid"                      "T_NK"                        
# [5] "B"                            "Basal_Myoepithelial"         
# [7] "Plasma"                       "Endothelial"                 
# [9] "CAF_adipocyte?"               "CAF_Perivascular"            
# [11] "Mast"                         "Epithelial (Lum Progenitor?)"
# [13] "pDC"                         
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT425B1-S1H1A3Y1N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                        "CAF_adipocyte?"              
# [3] "T_NK"                         "CAF"                         
# [5] "Myeloid"                      "Unknown_1"                   
# [7] "CAF_Perivascular"             "Plasma"                      
# [9] "Epithelial (Lum Progenitor?)" "Endothelial"                 
# [11] "B"                            "pDC"                         
# [13] "Basal_Myoepithelial"          "Mast"                        
# [15] "Epithelial (Tumor?)"         
# [1] "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/PM1467P1-T1Y2N1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv"
# [1] "Tumor"                  "Plasma"                 "Proliferating_T_cells" 
# [4] "CD4+T"                  "CD8+T"                  "Treg"                  
# [7] "TAM"                    "NK"                     "Pericyte"              
# [10] "myCAF"                  "apCAF"                  "mregDC"                
# [13] "AXL_DC"                 "B_cell"                 "iCAF"                  
# [16] "Endothelial"            "Monocyte"               "preAP_mesothelial_CAF" 
# [19] "CXCR4+CAF"              "Mast"                   "Adipocyte"             
# [22] "pDC"                    "cDC2"                   "cDC1"                  
# [25] "ILC"                    "Proliferating_myeloid"  "Vascular_smooth_muscle"
# HT065B1-S1H1A2A1Z1B1 and HT065B1-S1H7A2A1Z1B1 are not fully cell typed and have NA values in each of the cell type columns. 
# Additionally HT243B1-H3A2Fc1N1Z1_1Bn1_1 does not have coverage of the BRCA2 mutation site at chr13:32337161-32337164
# Because there are other sn/scRNA samples for these cases we will only proceed with the following single cell samples
# ccRCC matching object
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
obj <- readRDS("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/output/RCC_make_seurat_v1/CPT0079400003/RCC_CPT0079400003_RNA_and_SCT_0.3_resolution.rds")
sample.id <- "CPT0079400003"
cell_type_key <- "celltype_final"
obj$sample.id <- sample.id
# object has already been doublet filtered
obj$predicted_doublet <- "False" #every value in this column was "TRUE". Ilya confirmed that this is actually how he marks singlets. Relabeling as "False" to be consistent with the rest of the cohort data set.
obj <- UpdateSeuratObject(obj)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
meta.data <- read.table("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/RCC_clinical_v11_with_celltypes1.tsv", sep="\t", header=T)
sample_cell_type_meta <- meta.data[meta.data$orig.ident == sample.id, ]
rownames(sample_cell_type_meta) <- sample_cell_type_meta$merged_barcode
# tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = sample_cell_type_meta$merged_barcode)
# head(tmp_df)
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
#obj <- AddMetaData(obj, tmp_df , col.name = "cell_type_xenium_variant_v1")
obj <- AddMetaData(obj, sample_cell_type_meta)
obj$cell_type_xenium_variant_v1 <- obj$celltype_final
obj$celltype_final <- NULL
obj$celltype_final_short <- NULL
obj <- subset(x = obj, subset = cell_type_xenium_variant_v1 != 'doublets REMOVE')
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
all.genes <- rownames(obj)
print("starting normalization")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
print("plotting world domination")
pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
# print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("ccRCC cancer cell"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("ccRCC cancer cell")))] <- "normal" # there are no unknown cell labels in this object
saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
# conda activate r4-sig
library(Seurat)
library(Signac)
library(tidyverse)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/"
obj <- readRDS("/diskmnt/Projects/GBM_sc_analysis/multiome_autoprocess/autoprocess_objs/C3N-00663_CPT0087730015_2021-03-08/C3N-00663_CPT0087730015_2021-03-08_multiome_processed.rds")
obj <- UpdateSeuratObject(obj)
atac_codes <- colnames(obj@assays$ATAC@counts)
rna_codes <- colnames(obj@assays$RNA@counts)
all(atac_codes == rna_codes)
# [1] TRUE
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
saveRDS(obj,paste0(out_dir,"C3N-00663_CPT0087730015_2021-03-08_RNA_only_seuratv4.0.1.rds"))
# conda activate seurat5
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/C3N-00663_CPT0087730015_2021-03-08_RNA_only_seuratv4.0.1.rds")
sample.id <- "C3N-00663_CPT0087730015_2021-03-08"
obj <- UpdateSeuratObject(obj)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
cell_type_key <- "CellType"
obj$sample.id <- sample.id
meta.data <- read.table("/diskmnt/Projects/GBM_sc_analysis/scao/GBMLGG/1.integration/cptac.chop.gbm.lgg/2.xy/merged_40cptac_10chop_6lgg.rds.formated.cohort.ct.refined.121622.v2.tsv", sep="\t", header=T)
sample_cell_type_meta <- meta.data[meta.data$sample == sample.id, ]
rownames(sample_cell_type_meta) <- sample_cell_type_meta$merged_barcode
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
sample_cell_type_meta$orig.ident <- NULL
rownames(sample_cell_type_meta) <- sample_cell_type_meta$Barcode
sample_cell_type_meta$Barcode <- NULL
obj <- AddMetaData(obj, sample_cell_type_meta)
obj$cell_type_xenium_variant_v1 <- obj$CellType
obj$CellType <- NULL
# removing doublets
doublet_table <- "/diskmnt/Projects/GBM_sc_analysis/scao/GBMLGG/5.doublets/combo/combo/C3N-00663_CPT0087730015_2021-03-08/C3N-00663_CPT0087730015_2021-03-08_combo_scrublet_output_table.csv"
scrublet.df <- read.table(file = doublet_table, sep=",", header = TRUE, row.names = "Barcodes")
scrublet.df$predicted_doublet <- "NA"
scrublet.df$predicted_doublet[scrublet.df$predicted_doublet_rna == 'True' & scrublet.df$predicted_doublet_atac == 'True'] <- 'True'
scrublet.df$predicted_doublet[scrublet.df$predicted_doublet_rna == 'False' | scrublet.df$predicted_doublet_atac == 'False'] <- 'False'
new_row_names <- paste0(sample.id,"_",rownames(scrublet.df))
rownames(scrublet.df) <- new_row_names
obj <- AddMetaData(object = obj, metadata = scrublet.df)
obj <- subset(x = obj, predicted_doublet == 'False'| is.na(obj$predicted_doublet))
# adding metadata for consistent normalization
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
# normalizing
all.genes <- rownames(obj)
print("starting normalization")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
print("plotting world domination")
pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_rna", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_atac", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_rna", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_atac", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal" # there are no unknown cell labels in this object
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Unknown"))] <- "Unknown"
saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
# conda activate r4-sig
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original
library(Seurat)
library(Signac)
library(tidyverse)
set.seed(1234)
obj <- readRDS("/diskmnt/Projects/GBM_sc_analysis/multiome_autoprocess/autoprocess_objs/C3L-03372_CPT0275960013_2022-01-31/C3L-03372_CPT0275960013_2022-01-31_multiome_processed.rds")
obj <- UpdateSeuratObject(obj)
atac_codes <- colnames(obj@assays$ATAC@counts)
rna_codes <- colnames(obj@assays$RNA@counts)
all(atac_codes == rna_codes)
# [1] TRUE
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
saveRDS(obj,"/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/C3L-03372_CPT0275960013_2022-01-31_RNA_only_seuratv4.0.1.rds")
# conda activate seurat5
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/C3L-03372_CPT0275960013_2022-01-31_RNA_only_seuratv4.0.1.rds")
sample.id <- "C3L-03372_CPT0275960013_2022-01-31"
obj <- UpdateSeuratObject(obj)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
cell_type_key <- "cell_type"
obj$sample.id <- sample.id
meta.data <- read.table("/diskmnt/Projects/GBM_sc_analysis/cliu/Brain_tumor_evolution/Result/1_sample_processing/combo/C3L-03372_CPT0275960013_2022-01-31/C3L-03372_CPT0275960013_2022-01-31_metadata.tsv", sep="\t", header=T)
sample_cell_type_meta <- meta.data
sample_cell_type_meta$original_RNA_barcode <- rownames(sample_cell_type_meta)
sample_cell_type_meta$merged_barcode <- paste0(sample.id,"_",sample_cell_type_meta$original_RNA_barcode)
rownames(sample_cell_type_meta) <- sample_cell_type_meta$merged_barcode
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
sample_cell_type_meta$orig.ident <- NULL
sample_cell_type_meta$merged_barcode <- NULL
obj <- AddMetaData(obj, sample_cell_type_meta)
obj$cell_type_xenium_variant_v1 <- obj$cell_type
obj$cell_type <- NULL
# removing doublets
doublet_table <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/combo/C3L-03372_CPT0275960013_2022-01-31/C3L-03372_CPT0275960013_2022-01-31_combo_scrublet_output_table.csv"
scrublet.df <- read.table(file = doublet_table, sep=",", header = TRUE, row.names = "Barcodes")
scrublet.df$predicted_doublet <- "NA"
scrublet.df$predicted_doublet[scrublet.df$predicted_doublet_rna == 'True' & scrublet.df$predicted_doublet_atac == 'True'] <- 'True'
scrublet.df$predicted_doublet[scrublet.df$predicted_doublet_rna == 'False' | scrublet.df$predicted_doublet_atac == 'False'] <- 'False'
new_row_names <- paste0(sample.id,"_",rownames(scrublet.df))
rownames(scrublet.df) <- new_row_names
obj <- AddMetaData(object = obj, metadata = scrublet.df)
obj <- subset(x = obj, predicted_doublet == 'False'| is.na(obj$predicted_doublet))
# adding metadata for consistent normalization
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
# normalizing
all.genes <- rownames(obj)
print("starting normalization")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
print("plotting world domination")
pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_rna", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_atac", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_rna", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_atac", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor")))] <- "normal" # there are no unknown cell labels in this object
saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

# Subsetting out the matching scRNAseq from the Spatial PDAC driver paper cohort. Object and annotations were provided by Reyka.
#!/usr/bin/env Rscript
# HT060P1-S1R1A1G1Z1_1B2_1 -  cell barcode listed as HT060P_S1R1
# HT061P1-S1PAA1A1Z1_1B2_1 -  cell barcode listed as HT061P_S1PA
# HT125P1-S1H4A2K1G1Z1_1B1 - cell barcode listed as HT125P1_S1H4
# HT125P1-S1H8A2K1G1Z1_1B1 - cell barcode listed as HT125P1_S1H8
# HT185P1-S1H2A2K1G1Z1_1B1 - cell barcode listed as HT185P1_S1H2
# HT185P1-S1H3A2K1G1Z1_1B1 - cell barcode listed as HT185P1_S1H3
# /diskmnt/Projects/HTAN_analysis_2/PDAC_Pilot_Additional_Space/Additional_5/Everything_merge_v6.rds		
# /diskmnt/Projects/HTAN_analysis/PDAC/Spatial_PDAC_Metadata_Celltypes.txt
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
sample.vec <- c("HT060P1-S1R1A1G1Z1_1B2_1","HT061P1-S1PAA1A1Z1_1B2_1","HT125P1-S1H4A2K1G1Z1_1B1","HT125P1-S1H8A2K1G1Z1_1B1","HT185P1-S1H2A2K1G1Z1_1B1","HT185P1-S1H3A2K1G1Z1_1B1")
orig.ident.vec <- c("TWCE-HT060P1-S1R1A1A1Z1B1","TWCE-HT061P1-S1PAA1A1Z1B1","TWCE-HT125P1-XB2","TWCE-HT125P1-XB3","TWCE-HT185P1-XBc1_1","TWCE-HT185P1-XBc1_2")
obj_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC_Pilot_Additional_Space/Additional_5/Everything_merge_v6.rds"
cell_type_input = "/diskmnt/Projects/HTAN_analysis/PDAC/Spatial_PDAC_Metadata_Celltypes.txt"
cell_type_key = "cell_type_published"
cell_type_meta <- read.table(cell_type_input, sep='\t', header=T)
colnames(cell_type_meta) <- c("piece_barcode", cell_type_key)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
obj <- UpdateSeuratObject(obj)
for (i in 1:length(sample.vec)) {
    sample = sample.vec[i]
    orig.ident.sample = orig.ident.vec[i]
    print(sample)
    scsub <- subset(obj, subset = orig.ident == orig.ident.sample)
    
    sample_cell_type_meta <- cell_type_meta
    scsub$original_RNA_barcode <- stringr::str_sub(rownames(scsub@meta.data), -18, -1)
    tmp_df <- data.frame(cell_type_xenium_variant_v1 = sample_cell_type_meta[, cell_type_key], row.names = sample_cell_type_meta$piece_barcode)
    head(tmp_df)
    scsub <- AddMetaData(scsub, tmp_df , col.name = "cell_type_xenium_variant_v1")
    
    scrublet.df <- read.table(paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241111/RNA/",sample,"/",sample,"_scrublet_output_table.csv"), sep=",", header = TRUE, row.names = "Barcodes")
    scsub <- RenameCells(scsub, new.names = paste0(sample,"_",scsub$original_RNA_barcode))
    rownames(scrublet.df) <- paste0(sample,"_",rownames(scrublet.df))
    scsub <- AddMetaData(object = scsub, metadata = scrublet.df)
    
    scsub <- subset(x = scsub, subset = predicted_doublet == 'False'| is.na(scsub$predicted_doublet) | cell_type_xenium_variant_v1 == "ADM")
    scsub$percent.mito <- PercentageFeatureSet(scsub, pattern = "^MT-") / 100
    scsub$percent.rb <- PercentageFeatureSet(scsub, pattern = "^RBS|RPL", assay = "RNA")
    all.genes <- rownames(scsub)
    scsub <- NormalizeData(scsub, normalization.method = "LogNormalize", scale.factor = 10000)
    scsub <- FindVariableFeatures(scsub, selection.method = "vst", nfeatures = 2000)
    scsub <- ScaleData(scsub, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
    scsub <- SCTransform(scsub, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(scsub)), return.only.var.genes = FALSE)
    DefaultAssay(scsub) <- "SCT"
    scsub <- RunPCA(scsub, npcs = 200, verbose = TRUE)
    scsub <- RunUMAP(scsub, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
    scsub <- FindNeighbors(scsub, reduction = "pca", dims = 1:30, force.recalc = T)
    scsub <- FindClusters(scsub, resolution = 0.5)
    Idents(scsub) <- "seurat_clusters"
    pdf(paste0(out_dir,"/","UMAP_QC_",sample,".pdf"), useDingbats = F)
    print(rasterize(DimPlot(scsub, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(scsub, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(scsub, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(scsub, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(scsub, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
    print(rasterize(FeaturePlot(scsub, reduction = "umap.30PC", features = "doublet_score", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(scsub, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(scsub, reduction = "umap.30PC", group.by = "KRAS_hotspots2", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    print(rasterize(DimPlot(scsub, reduction = "umap.30PC", group.by = "SMGs2", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    pdf(paste0(out_dir,"/","UMAP_cell_type_xenium_varaint_v1",sample,"_clara.pdf"))
    print(rasterize(DimPlot(scsub, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
    dev.off()
    print(table(scsub$cell_type_xenium_variant_v1, scsub$seurat_clusters))
    scsub$neoplasm_normal_unknown <- NA
    scsub$neoplasm_normal_unknown[(scsub$cell_type_xenium_variant_v1 %in% c("PDAC","PanIN"))] <- "neoplasm"
    scsub$neoplasm_normal_unknown[(!(scsub$cell_type_xenium_variant_v1 %in% c("PDAC","PanIN")))] <- "normal"
    saveRDS(scsub, paste0(out_dir,"/",sample,"_processed_update_no_doublet_seurat5.0.1.rds"))
    tmp_df <- scsub@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
    write.table(tmp_df, paste0(out_dir,"/",sample,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    tmp_df <- scsub@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
    write.table(tmp_df, paste0(out_dir,"/",sample,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}
# conda activate r4-sig
library(Seurat)
library(Signac)
library(tidyverse)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/"
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/TNP_CRC/seurat/combo_v5.0/HT179C1-T1A3Y1N1/HT179C1-T1A3Y1N1_processed_multiomic.rds")
obj <- UpdateSeuratObject(obj)
atac_codes <- colnames(obj@assays$ATAC@counts)
rna_codes <- colnames(obj@assays$RNA@counts)
all(atac_codes == rna_codes)
# [1] TRUE  
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
saveRDS(obj,paste0(out_dir,"HT179C1-T1A3Y1N1_RNA_only_seuratv4.0.1.rds")) #does not have doublets removed so need to run scrublet as well.
# checking markers for the cell typed single cell objects
# conda activate seurat5
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/"
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/3_Merged_scRNA_179C1_TNP_CRC_092121.rds")
sample.id <- "3_Merged_scRNA_179C1_TNP_CRC_092121"
obj <- UpdateSeuratObject(obj)
# cell_type_key <- "CellType"
# obj$sample.id <- sample.id
meta.data <- read.table("/diskmnt/Projects/HTAN_analysis_2/TNP_CRC/cell_type/v5.0/manual/179C1_scRNA/179C1_scRNA_processed_multiomic_cellTyped.meta.data", sep="\t", header=T)
rownames(meta.data) <- meta.data$X
tmp <- data.frame(cell_types = meta.data$cell_type, row.names = rownames(meta.data))
obj <- AddMetaData(obj, tmp, col.name = "cell_type")
unique(obj$cell_type)
library(future)
options(future.globals.maxSize= +Inf)
plan(multicore, workers = 50)
Idents(obj) <- "cell_type"
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
write.table(degs,paste0(out_dir,"/",sample.id,"_FindAllMarkers_cell_type_20241130.tsv"),sep="\t",quote=FALSE)
plan(sequential)
obj$barcode <- rownames(obj@meta.data)
non_doublet_cells <- obj$barcode[(!(obj$cell_type == "Doublets"))]
sobj <- subset(x = obj, subset = barcode %in% c(non_doublet_cells))
sobj$percent.rb <- PercentageFeatureSet(sobj, pattern = "^RBS|RPL", assay = "RNA")
all.genes <- rownames(sobj)
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
sobj <- ScaleData(sobj, features = all.genes, vars.to.regress = c('percent.mt', 'percent.rb','S.Score','G2M.Score'))
sobj <- SCTransform(sobj, vst.flavor="v2", vars.to.regress = c('percent.mt', 'percent.rb','S.Score','G2M.Score'), ncells = length(Cells(sobj)), return.only.var.genes = FALSE)
DefaultAssay(sobj) <- "SCT"
sobj <- RunPCA(sobj, npcs = 200, verbose = TRUE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:30, force.recalc = T)
sobj <- FindClusters(sobj, resolution = 0.5)
Idents(sobj) <- "cell_type"
pdf(paste0(out_dir,"UMAP_QC_Doublets_removed_",sample.id,".pdf"), useDingbats = F)
print(rasterize(DimPlot(sobj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(sobj, reduction = "umap.30PC", group.by = "cell_type", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(sobj, reduction = "umap.30PC", features = c("percent.mt"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(sobj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(sobj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(sobj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
dev.off()
plan(multicore, workers = 50)
Idents(sobj) <- "cell_type"
degs <- FindAllMarkers(sobj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
write.table(degs,paste0(out_dir,"/",sample.id,"_FindAllMarkers_cell_type_no_Doublets_20241130.tsv"),sep="\t",quote=FALSE)
# conda activate seurat5
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/original/HT179C1-T1A3Y1N1_RNA_only_seuratv4.0.1.rds")
sample.id <- "HT179C1-T1A3Y1N1"
obj <- UpdateSeuratObject(obj)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = c('RNA'))
cell_type_key <- "CellType"
obj$sample.id <- sample.id
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
obj$barcode <- paste0(sample.id,"_",obj$original_RNA_barcode)
obj <- RenameCells(obj, new.names = paste0(sample.id,"_",obj$original_RNA_barcode))
# removing doublets
doublet_table <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/scrublet/20241130/combo/HT179C1-T1A3Y1N1/HT179C1-T1A3Y1N1_combo_scrublet_output_table.csv"
scrublet.df <- read.table(file = doublet_table, sep=",", header = TRUE, row.names = "Barcodes")
scrublet.df$predicted_doublet <- "NA"
scrublet.df$predicted_doublet[scrublet.df$predicted_doublet_rna == 'True' & scrublet.df$predicted_doublet_atac == 'True'] <- 'True'
scrublet.df$predicted_doublet[scrublet.df$predicted_doublet_rna == 'False' | scrublet.df$predicted_doublet_atac == 'False'] <- 'False'
new_row_names <- paste0(sample.id,"_",rownames(scrublet.df))
rownames(scrublet.df) <- new_row_names
obj <- AddMetaData(object = obj, metadata = scrublet.df)
obj <- subset(x = obj, predicted_doublet == 'False'| is.na(obj$predicted_doublet))
# adding metadata for consistent normalization
obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-") / 100
obj$percent.rb <- PercentageFeatureSet(obj, pattern = "^RBS|RPL", assay = "RNA")
# normalizing
all.genes <- rownames(obj)
print("starting normalization")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes, vars.to.regress = c('percent.mito', 'percent.rb'))
obj <- SCTransform(obj, vst.flavor="v2", vars.to.regress = c('percent.mito', 'percent.rb'), ncells = length(Cells(obj)), return.only.var.genes = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 200, verbose = TRUE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.30PC", reduction.key = "UMAP30PC_")
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)
Idents(obj) <- "seurat_clusters"
print("plotting world domination")
pdf(paste0(out_dir,"UMAP_QC_",sample.id,".pdf"), useDingbats = F)
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "seurat_clusters", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.mito"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("percent.rb"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nFeature_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = c("nCount_RNA"), raster = F),layers='Point', dpi=300))
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_rna", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(FeaturePlot(obj, reduction = "umap.30PC", features = "doublet_score_atac", label=TRUE, label.size=4, raster=FALSE),layers='Point', dpi=300)) # this field is not in the metadata
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_rna", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet_atac", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "predicted_doublet", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
library(future)
options(future.globals.maxSize= +Inf)
plan(multicore, workers = 50)
Idents(obj) <- "seurat_clusters"
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
write.table(degs,paste0(out_dir,"/",sample.id,"_FindAllMarkers_seurat_clusters_20241130.tsv"),sep="\t",quote=FALSE)
plan(sequential)
saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/marker_dotplots/20241202
# Rscript $austin/tools/markers/cell_typing_plots_Liver.R .
# Rscript $austin/tools/markers/cell_typing_plots_PDAC_caf.R .
# Rscript $austin/tools/markers/cell_typing_plots.R .
# 0 Tumor
# 1 Tumor
# 2 Tumor
# 3 Tumor
# 4 Myeloid
# 5 Tumor_proliferative
# 6 iCAF
# 7 Hepatocytes 
# 8 T_cell
# 9 Endothelial
# 10 myCAF
# 11 Tumor center of tumor clusters marked by high mito reads
# 12 Tumor
# 13 Cholangiocytes
# 14 Plasma <- probably plasma based on IGHG genes expression but could also include B_cell and pDC
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1//HT179C1-T1A3Y1N1_processed_update_no_doublet_seurat5.0.1.rds
obj$cell_type_xenium_variant_v1 <- NA
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(0,1,2,3,11,12)] = "Tumor"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(4)] = "Myeloid"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(5)] = "Tumor_proliferative"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(6)] = "iCAF"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(7)] = "Hepatocytes"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(8)] = "T_cell"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(9)] = "Endothelial"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(10)] = "myCAF"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(13)] = "Cholangiocytes"
obj$cell_type_xenium_variant_v1[obj$seurat_clusters %in% c(14)] = "Plasma"
pdf(paste0(out_dir,"UMAP_cell_type_xenium_varaint_v1",sample.id,".pdf"))
print(rasterize(DimPlot(obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1", label=TRUE, label.size=2, raster=FALSE),layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_xenium_variant_v1, obj$seurat_clusters))
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type_xenium_variant_v1 %in% c("Tumor","Tumor_proliferative"))] <- "neoplasm"
obj$neoplasm_normal_unknown[(!(obj$cell_type_xenium_variant_v1 %in% c("Tumor","Tumor_proliferative")))] <- "normal" # there are no unknown cell labels in this object
saveRDS(obj, paste0(out_dir,"/",sample.id,"_processed_update_no_doublet_seurat5.0.1.rds"))
tmp_df <- obj@meta.data[,c("original_RNA_barcode","cell_type_xenium_variant_v1")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_cell_type_xenium_variant_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write.table(tmp_df, paste0(out_dir,"/",sample.id,"_neoplasm_normal_unknown_v1.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)