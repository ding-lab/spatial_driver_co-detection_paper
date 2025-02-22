#!/usr/bin/env Rscript
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/
# Rscript plotting_neoplastic_labels_7_darkblue_ceiling1.R neo_norm_unk_table_v7.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table_path <- args[1]
out_dir <- args[2]
#out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/"
#input_table <- read.table("neo_norm_unk_table_v3.tsv",sep='\t',header=T,row.names="X")
input_table <- read.table(input_table_path,sep='\t',header=T)
# Sample_ID	Xenium_neoplastic_labels_v6	Xenium_neoplastic_labels_v6_subclone	Xenium_unknown_low_quality_labels_v6	Manual_xenium_cell_types_v6	Xenium_snv_object	feature_csv	mutant_csv	feature_csv_table_plots
iterables <- rownames(input_table)
for (i in iterables) {
    sample = input_table[i,"Sample_ID"]
    rds_obj = input_table[i,"Xenium_snv_object"]
    dir.create(paste0(out_dir,"/",sample))
    Manual_cell_type = input_table[i,"Manual_xenium_cell_types_v6"]
    neoplastic_labels_v3 <- input_table[i,"Xenium_neoplastic_labels_v6"]
    neoplastic_labels_v3_subclone <- input_table[i,"Xenium_neoplastic_labels_v6_subclone"]
    unknown_low_quality_labels_v3 <- input_table[i,"Xenium_unknown_low_quality_labels_v6"]
    tumor_labels <- unlist(strsplit(neoplastic_labels_v3, ","))
    tumor_labels_subclone <- unlist(strsplit(neoplastic_labels_v3_subclone, ","))
    unknown_low_quality_labels_v3 <- unlist(strsplit(unknown_low_quality_labels_v3, ","))
    feature_csv = input_table[i,"feature_csv"]
    feature_list <- unlist(strsplit(feature_csv, ","))
    print(sample)
    print(feature_list)
    obj <- readRDS(rds_obj)
    DefaultAssay(obj) <- "Xenium.with.snvs"
    Idents(obj) <- "seurat_clusters"
    DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
    DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
    cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
    obj$barcode <- rownames(obj@meta.data)
    colnames(cell_types) <- c("barcode","cell_type")
    rownames(cell_types) <- cell_types$barcode
    cell_types$barcode <- NULL
    obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
    obj$neoplasm_normal_unknown <- NA
    obj$neoplasm_normal_unknown[(obj$cell_type %in% tumor_labels)] <- "neoplastic"
    obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "normal"
    obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
    # neoplasm_normal_unknown color codes c('#CC79A7','#009E73','#B4DAF4') ordered by normal, neoplasm, unknown
    obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic","low_quality/unknown"))
    tmp_df <- obj@meta.data[,c("barcode","neoplasm_normal_unknown")]
    colnames(tmp_df) <- c("cell_id","neoplasm_normal_unknown_v3")
    write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplasm_normal_unknown_v4.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
    tmp_df <- obj@meta.data[,c("barcode","neoplasm_normal_unknown")]
    colnames(tmp_df) <- c("cell_id","group")
    write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplasm_normal_unknown_v4.csv"), row.names=FALSE, sep=",", quote=FALSE)
    obj$uniform_background_color = "darkblue"
    # the featureplot function first finds the xlims and ylims and then uses them to set the limits and then can do coord_fixed (by default coord_fixed is set to false) and then does some stuff regarding the legend scaling and if there is more than one featureplot and then returns that featureplot.
    # xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
    # ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
    # ...
    # +
    #     scale_x_continuous(limits = xlims) +
    #     scale_y_continuous(limits = ylims) +
    # ...
    # # Add coord_fixed
    # if (coord.fixed) {
    #     plot <- plot + coord_fixed()
    # }
    # # I'm not sure why, but sometimes the damn thing fails without this
    # # Thanks ggplot2
    # plot <- plot
    # 
    coordinates <- Embeddings(obj, reduction = "umap.30PC")
    dim_names <- colnames(coordinates)
    dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
    dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
    sorted_dim_ylims <- sort(dim_ylims)
    sorted_dim_xlims <- sort(dim_xlims)
    x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
    y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
    xy_ratio <- x_length/y_length
    #variant_transcript = "#59B3E4"
    #reference_transcript = "#EFE642"
    #from the Okabe and Ito color palette https://siegal.bio.nyu.edu/color-palette/ https://jfly.uni-koeln.de/color/
    new_variant_transcript = "#56B4E9"
    new_reference_transcript = "F0E442"
    p1 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), reduction = "umap.30PC", pt.size=2.5, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
    p2 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=2.5, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
    list_of_plots <- list(p1)
    jiterables = length(list_of_plots)
    list_of_plots[[ length(list_of_plots) + 1 ]] <- p2
    jiterables = length(feature_list)
    for (j in 1:jiterables) {
        p3 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], max.cutoff = 1, order = T, reduction = "umap.30PC", label=F, pt.size=2.5, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
        list_of_plots[[j+2]] <- p3
    }
    print("Feature plot length")
    print(length(list_of_plots))
    pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4.pdf"),useDingbats = F,height=10, width = 10*(length(list_of_plots)))
    print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
    dev.off()
    # pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
    # for (j in 1:length(list_of_plots)) {
    #     print(list_of_plots[[j]])
    # }
    # dev.off()
    p4 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), reduction = "umap.30PC", pt.size=2.5, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
    p5 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=2.5, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + guides(color=guide_legend(ncol=1, override.aes = list(size=3))), layers="Point", dpi=300)
    list_of_plots <- list(p4)
    list_of_plots[[ length(list_of_plots) + 1 ]] <- p5
    jiterables = length(feature_list)
    for (j in 1:jiterables) {
        p6 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], max.cutoff = 1, order = T, reduction = "umap.30PC", label=F, pt.size=2.5, raster=FALSE) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
        list_of_plots[[j+2]] <- p6
    }
    print("Feature plot legend length")
    print(length(list_of_plots))
    pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_legend.pdf"),useDingbats = F,height=10, width=10*(length(list_of_plots)))
    print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
    dev.off()
    # pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
    # for (j in 1:length(list_of_plots)) {
    #     print(list_of_plots[[j]])
    # }
    # dev.off()
    # Image 3 color without legend
    #p7 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    #p8 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    #list_of_plots <- list(p7)
    #list_of_plots[[ length(list_of_plots) + 1 ]] <- p8
    #jiterables = length(feature_list)
    #for (j in 1:jiterables) {
    #    p9 <- rasterize(ImageFeaturePlot(obj, features = feature_list[j], max.cutoff = 'q90', size = 0.15, cols = c("#141414", "darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    #    list_of_plots[[j+2]] <- p9
    #}
    #print("Image Feature plot 3 length")
    #print(length(list_of_plots))
    #pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_3color.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
    #print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
    #dev.off()
    # image 3 color with legend
    #p10 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    #p11 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    #list_of_plots <- list(p10)
    #list_of_plots[[ length(list_of_plots) + 1 ]] <- p11
    #jiterables = length(feature_list)
    #for (j in 1:jiterables) {
    #    p12 <- rasterize(ImageFeaturePlot(obj, features = feature_list[j], max.cutoff = 'q90', size = 0.15, cols = c("#141414", "darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    #    list_of_plots[[j+2]] <- p12
    #}
    #print("Image Feature plot 3 length legend")
    #print(length(list_of_plots))
    #pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_3color_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
    #print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
    #dev.off()
    # Image 3 color without legend
    p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    list_of_plots <- list(p13)
    list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
    jiterables = length(feature_list)
    for (j in 1:jiterables) {
        p15 <- rasterize(ImageFeaturePlot(obj, features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
        list_of_plots[[j+2]] <- p15
    }
    p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    offset = length(list_of_plots)
    list_of_plots[[1+offset]] <- p15
    print("Image Feature plot 2 length")
    print(length(list_of_plots))
    pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_image_2color.pdf"),useDingbats = F,height=10, width = 10*(length(list_of_plots)))
    print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
    dev.off()
    # image 2 color with legend
    p16 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    p17 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) + guides(color=guide_legend(ncol=1, override.aes = list(size=4L, alpha=1))), layers="Polygon", dpi=600)
    list_of_plots <- list(p16)
    list_of_plots[[ length(list_of_plots) + 1 ]] <- p17
    jiterables = length(feature_list)
    for (j in 1:jiterables) {
        p18 <- rasterize(ImageFeaturePlot(obj, features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
        list_of_plots[[j+2]] <- p18
    }
    p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
    offset = length(list_of_plots)
    list_of_plots[[1+offset]] <- p15
    print("Image Feature plot 2 length legend")
    print(length(list_of_plots))
    pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_image_2color_legend.pdf"),useDingbats = F,height=10, width=10*(length(list_of_plots)))
    print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
    dev.off()
    # print(ImageFeaturePlot(obj, features = feature_list[j], max.cutoff = 'q95', size = 0.15, cols = c("#3b3b3b", "darkblue", "yellow"), border.color= NA)) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
    # print(ImageFeaturePlot(obj, features = feature_list[j], max.cutoff = 'q90', size = 0.15, cols = c("darkblue", "yellow"), border.color= NA)) #max.cutoff is 90th percentile
    #
    # pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
    # for (j in 1:length(list_of_plots)) {
    #     print(list_of_plots[[j]])
    # }
    # dev.off()
}
