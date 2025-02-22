#!/usr/bin/env Rscript
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/"
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table <- read.table(args[1],sep='\t',header=T,row.names="X")
out_dir = args[2]
# input_table <- read.table("plotting_input_table.tsv",sep='\t',header=T,row.names="X")
# sampleID, rds_obj, Manual_cell_type, variant, reference
iterables <- rownames(input_table)
for (i in iterables) {
    sample = input_table[i,"sampleID"]
    rds_obj = input_table[i,"rds_obj"]
    Manual_cell_type = input_table[i,"Manual_cell_type"]
    variant_table = input_table[i,"variant"]
    reference_table = input_table[i,"reference"]
    print(sample)
    obj <- readRDS(rds_obj)
    DefaultAssay(obj) <- "Xenium.with.snvs"
    Idents(obj) <- "seurat_clusters"
    DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
    DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
    variant_alleles <- read.table(variant_table,sep='\t',header=F)
    reference_alleles <- read.table(reference_table,sep='\t',header=F)
    cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
    colnames(cell_types) <- c("barcode","cell_type")
    rownames(cell_types) <- cell_types$barcode
    cell_types$barcode <- NULL
    obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
    pdf(paste0(out_dir,"/",sample,"/",sample,"_manual_cell_types_v1.pdf"),useDingbats = F, height=10, width = 14)
    print(DimPlot(obj, group.by = "seurat_clusters", reduction = "umap.30PC", label=T, label.size=3, raster=FALSE) + coord_fixed(ratio=1))
    print(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", label=T, label.size=3, raster=FALSE) + coord_fixed(ratio=1))
    dev.off()
    pdf(paste0(out_dir,"/",sample,"/",sample,"_genotype_featureplots.pdf"),useDingbats = F,height=5, width = 7)
    for (j in 1:length(variant_alleles[,"V1"])) {
        variant_key <- variant_alleles[j,'V1']
        wt_key <- reference_alleles[j,'V1']
        print(variant_key)
        print(wt_key)
        p1 <- FeaturePlot(obj, slot="counts", features = variant_key, order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=1)#, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
        p2 <- FeaturePlot(obj, slot="counts", features = wt_key, order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=1)#, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
        print(rasterize(p1, layers='Point', dpi=300))
        print(rasterize(p2, layers='Point', dpi=300))
    }
    dev.off()
}
