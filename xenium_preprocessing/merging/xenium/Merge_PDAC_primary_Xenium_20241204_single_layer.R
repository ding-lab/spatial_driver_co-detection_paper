##!/usr/bin/env Rscript
##################################################################
## Merge Xenium objects with variant probes and multiple panels ##
## Author: Austin Southard-Smith                                ##
## Last edited: 10/23/2023                                      ##
## Environment: conda activate seurat5                          ##
##################################################################
# conda activate seurat5
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(Seurat)
library(optparse)
library(patchwork)
set.seed(1234)
library(future)
options(future.globals.maxSize= +Inf)
#plan(multicore, workers = 90)
plan(sequential)
# Rscript Merge_script.R -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/merge_input_test_table.tsv -p PDAC_merge_primary_KRAS_test -o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/test/
option_list = list(
    make_option(c("-i", "--input_table"),
                type="character",
                default=NULL,
                help="path to input_table tsv file containing paths to xenium objects,  with seurat objects to merge; each subdirectory should be as follows: sample_id/sample_id_processed.rds",
                metavar="character"),
    make_option(c("-o", "--output_dir"),
                type="character",
                default=NULL,
                help="path to desired output directory where merged objects will be written to",
                metavar="character"),
    make_option(c("-p", "--prefix"),
                type="character",
                default=NULL,
                help="prefix for output files (e.g. snRNA_BRCA, BRCA, HTAN_BRCA)",
                metavar="character"),
    make_option(c("--counts_cutoff"),
                type="integer",
                default=0,
                help="retain only cells with nCounts greater than the value specified by the counts_cutoff",
                metavar="integer"),
    make_option(c("--pc_num_individual"),
                type="integer",
                default=30,
                help="number of principal components to use for individual PCA calculation post feature subsetting",
                metavar="integer")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt$input_table <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/merge_input_test_table.tsv"
#opt$prefix <- "PDAC_merge_primary_KRAS_test"
#opt$output_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/test/"
if (is.null(opt$input_table)){
    print_help(opt_parser)
    stop("Path to data is required (-i; --input_table).n", call.=FALSE)
}
if (is.null(opt$output_dir)){
    print_help(opt_parser)
    stop("Path to output directory is required (-o; --output_dir).n", call.=FALSE)
}
if (is.null(opt$prefix)){
    print_help(opt_parser)
    stop("Prefix for output files is required (-r; --prefix).n", call.=FALSE)
}

out_dir = opt$output_dir
dir.create(out_dir)
setwd(out_dir)
prefix = opt$prefix
# all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v5_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
# input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
input_table <- read.table(opt$input_table, sep='\t',header=T)
# out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1/"
# sample_ID piece_ID    case_ID obj_path    Manual_xenium_cell_types_v4  Xenium_neoplastic_labels_v6 Xenium_neoplastic_labels_v6_subclone    Xenium_unknown_low_quality_labels_v6    Progression Sex Site    Disease   Treatment   Age   Section_KRAS_call
pc_num_individual = opt$pc_num_individual
print("Loading objects to merge ...")
input_table <- read.table(opt$input_table, sep='\t',header=T)
seurat.list = c()
iterables <- 1:length(rownames(input_table))
sample_vector = c()
for (i in iterables) {
    obj_path = input_table[i,"Xenium_snv_object"]
    sampleID = input_table[i,"sample_ID"]
    sample_vector = c(sample_vector, sampleID)
    print(sampleID)
    caseID = input_table[i,"case_ID"]
    pieceID = input_table[i,"piece_ID"]
    progression = input_table[i,"Progression"]
    sex = input_table[i,"Sex"]
    site = input_table[i,"Site"]
    disease = input_table[i,"Disease"]
    treatment = input_table[i,"Treatment"]
    age = input_table[i,"Age"]
    section_KRAS_call = input_table[i,"Section_KRAS_call"]
    neoplastic_labels_v4 <- input_table[i,"Xenium_neoplastic_labels_v6"]
    neoplastic_labels_v4_subclone <- input_table[i,"Xenium_neoplastic_labels_v6_subclone"]
    unknown_low_quality_labels_v4 <- input_table[i,"Xenium_unknown_low_quality_labels_v6"]
    tumor_labels <- unlist(strsplit(neoplastic_labels_v4, ","))
    tumor_labels_subclone <- unlist(strsplit(neoplastic_labels_v4_subclone, ","))
    unknown_low_quality_labels_v3 <- unlist(strsplit(unknown_low_quality_labels_v4, ","))
    obj_new = readRDS(obj_path)
    obj_new$sample_ID = sampleID
    obj_new$case_ID = caseID
    obj_new$piece_ID = pieceID
    obj_new$Progression = progression
    obj_new$Sex = sex
    obj_new$Site = site
    obj_new$Disease = disease
    obj_new$Treatment = treatment
    obj_new$Age = age
    obj_new$Section_KRAS_call = section_KRAS_call
    manual_cell_type_path = input_table[i,"Manual_xenium_cell_types_v6"]
    cell_types <- read.table(manual_cell_type_path, header = T, sep = '\t')
    obj_new$original_Xenium_barcode <- rownames(obj_new@meta.data)
    colnames(cell_types) <- c("barcode","cell_type")
    rownames(cell_types) <- cell_types$barcode
    cell_types$barcode <- NULL
    # Features(obj_new, assay = "Xenium.with.snvs")
    assay_list <- Assays(obj_new)
    assays_with_snv_probes <- c()
    for (assay in assay_list) {
        DefaultAssay(obj_new) <- assay
        feature_vector <- dimnames(obj_new)[1][[1]]
        #print("feature_vector")
        #print(feature_vector)
        if (sum((grepl("-ALT",feature_vector)) | (grepl("-WT", feature_vector))) > 0) {
            assays_with_snv_probes <- c(assays_with_snv_probes, assay)
        }
    }
    # transferring the snv counts to the metadata field because it will make working with integrated objects easier as we need to drop all non-shared genes from the object and the easiest way to do that is by subsetting.
    for (assay in assays_with_snv_probes) {
        print(assay)
        assay_class = class(obj_new[[assay]])
        if (assay_class == "Assay5") {
            #print(1)
            counts_df <- t(as.matrix(GetAssayData(object = obj_new, assay = assay, layer = "counts")))
            counts_df <- as.data.frame(counts_df)
            #print(counts_df[1:5,1:5])
        } else {
            counts_df <- as.data.frame(t(as.matrix(GetAssayData(object = obj_new, assay = assay, slot = "counts"))))
        }
        assay_feature_list <- Features(obj_new, assay=assay)
        variant_key_list <- assay_feature_list[grepl("-ALT",assay_feature_list)]
        wt_key_list <- assay_feature_list[grepl("-WT",assay_feature_list)]
        num_var_keys <- length(variant_key_list)
        num_ref_keys <- length(wt_key_list)
        for (j in 1:num_var_keys) {
            variant_key <- variant_key_list[j]
            print(variant_key)
            var_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
            # making shared variant names consistent
            if (variant_key == "PIK3CA-p-E545K-ALT-21-A") { # This gene is on the positive strand so it is simila
                variant_key <- "PIK3CA-p-E545K-ALT-A"
            } else if (variant_key == "KRAS-p-G12V-ALT-21-T") { # The allele here of T is based on the opposite strand (gotta love negative vs positive naming conventions) of the one from the new name. The actual probe sequence used is identical.
                variant_key <- "KRAS-p-G12V-ALT-A"
            } else if (variant_key == "KRAS-p-G12D-ALT-21-A") { # The allele here of A is based on the opposite strand (gotta love negative vs positive naming conventions) of the one from the new name. The actual probe sequence used is identical.
                variant_key <- "KRAS-p-G12D-ALT-T"
            }
            obj_new <- AddMetaData(object = obj_new, metadata = var_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
        }
        for (j in 1:num_ref_keys) {
            wt_key <- wt_key_list[j]
            print(wt_key)
            wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
            obj_new <- AddMetaData(object = obj_new, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
        }
    }
    obj_new <- AddMetaData(object = obj_new, metadata = cell_types, col.name = "cell_type_individual")
    # obj_new <- RenameCells(obj_new, new.names = paste0(sampleID,"_",obj_new$original_Xenium_barcode)) # do this for integration rather than merging.
    obj_new$neoplasm_normal_unknown <- NA
    obj_new$neoplasm_normal_unknown[(obj_new$cell_type_individual %in% tumor_labels)] <- "neoplastic"
    obj_new$neoplasm_normal_unknown[(!(obj_new$cell_type_individual %in% tumor_labels))] <- "normal"
    obj_new$neoplasm_normal_unknown[(obj_new$cell_type_individual %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
    obj_new$neoplasm_normal_unknown_subclone <- NA
    obj_new$neoplasm_normal_unknown_subclone[(obj_new$cell_type_individual %in% tumor_labels_subclone)] <- "neoplastic"
    obj_new$neoplasm_normal_unknown_subclone[(!(obj_new$cell_type_individual %in% tumor_labels_subclone))] <- "normal"
    obj_new$neoplasm_normal_unknown_subclone[(obj_new$cell_type_individual %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
    obj_new$seurat_clusters_individual <- obj_new$seurat_clusters
    reduction_list <- Reductions(obj_new)
    for (reduc in reduction_list) {
        if (sum(grepl("UMAP", toupper(reduc)))>0) {
            coordinates <- Embeddings(obj_new, reduction = reduc)
            dim_names_new <- paste0(colnames(coordinates),"_individual")
            obj_new <- AddMetaData(obj_new, metadata = coordinates, col.name = dim_names_new)
        }
    }
    seurat.list = c(seurat.list, obj_new)
}

# Can't just do normal merging where features present in one object but not another are auto filled in with zeros because 
# unlike in single cell data if a gene is present in the matrix of one xenium seurat object but not another it is because that gene was not included on the panel and therefore not included for measurement.
# as a result we can only look at shared features here.
# finding the intersection of features from each Xenium object.
print("identifying shared feature set across all samples to filter to prior to merge/integration. Retaining all controls.")
intersecting_features <- c()
control_features <- c()
snvs_to_preserve <- c()
for (i in iterables) {
    obj_new <- seurat.list[[i]]
    sampleID <- unique(obj_new$sample_ID)
    obj_features <- Features(obj_new, assay = c("Xenium"))
    if (i == iterables[1]) {
        intersecting_features <- c(intersecting_features, obj_features)
    } else {
        intersecting_features <- intersect(intersecting_features, obj_features)
    }
    snv_features <- Features(obj_new, assay = "Xenium.with.snvs")
    snv_features <- snv_features[(grepl("-ALT",snv_features)) | (grepl("-WT", snv_features))]
    snvs_to_preserve <- c(snvs_to_preserve, snv_features)
    for (assay in c("BlankCodeword","ControlCodeword","ControlProbe")) {
        obj_controls <- Features(obj_new, assay = assay)
        control_features <- unique(c(control_features, obj_controls))
    }
}
# this is the same subset_opt function used to initially generate the objects. This file is present at `spatial_driver_co-detection_paper/xenium_preprocessing/seurat_object_generation/subset_obj_seurat.R` on github
source("/diskmnt/Projects/Users/austins2/tools/subset_obj_seurat.R")
print("Subsetting each individual sample to only a shared feature set and then re-normalizing each individual object")
total_cells = 0
for (i in iterables) {
    obj_new <- seurat.list[[i]]
    DefaultAssay(obj_new) <- "Xenium"
    sampleID = unique(obj_new$sample_ID)
    print(sampleID)
    obj_new$nCount_Xenium_individual <- obj_new$nCount_Xenium
    obj_new$nFeature_Xenium_individual <- obj_new$nFeature_Xenium
    obj_new$nCount_Xenium.with.snvs_individual <- obj_new$nCount_Xenium.with.snvs
    obj_new$nFeature_Xenium.with.snvs_individual <- obj_new$nFeature_Xenium.with.snvs
    obj_new$nCount_SCT_individual <- obj_new$nCount_SCT
    obj_new$nFeature_SCT_individual <- obj_new$nFeature_SCT
    # this option also updates the nCount_Xenium, nFeature_Xenium, nCount_Xenium.with.snvs, nFeature_Xenium.with.snvs, nCount_SCT, nFeature_SCT
    # we do not drop any control probe features so the following are not updated and need to be copied to a different location in the metadata: nCount_BlankCodeword, nFeature_BlankCodeword, nCount_ControlCodeword, nFeature_ControlCodeword, nCount_ControlProbe, nFeature_ControlProbe
    obj_new <- subset_opt(obj_new, features = c(intersecting_features, control_features, snvs_to_preserve))
    # because we are subsetting features some cells may no longer have any counts in them (i.e. nCounts_Xenium = 0) so we remove them based on the filtering conditions previously established for the study
    if (sum(obj_new$nCount_Xenium <= opt$counts_cutoff) > 0) {
        obj_new.empty <- subset_opt(obj_new, subset = nCount_Xenium <= opt$counts_cutoff)
        empty_cells = Cells(obj_new.empty)
    } else {
        empty_cells = c()
    }
    write.table(empty_cells, paste(out_dir,"/","Empty_cells_post_feature_intersection_subset_",sampleID,'.tsv',sep=''),sep="\t",quote=FALSE)
    rm(obj_new.empty)
    # this step removes the Xenium.with.snvs assay to save memory and storage space since it is now identical to the "Xenium" operation after the feature subset operation
    # right now it is not necessary but it may be necessary for work with larger integrated/merged datasets.
    # DietSeurat(obj_new, assays = c("BlankCodeword","ControlCodeword","ControlProbe","Xenium"))
    # DefaultAssay(obj_new) <- "Xenium"
    obj_new <- subset_opt(obj_new, subset = nCount_Xenium > opt$counts_cutoff)
    non_empty_cells <- Cells(obj_new)
    print("Filtering the image data to only the cells remaining in the Xenium object")
    image_names <- Images(obj_new)
    for (image_name in image_names) {
        print(paste0("updating FOV ",image_name))
        image_new <- obj_new[[image_name]]
        # the above subset operation already removes cells from the counts matrices and the fovs so no need to remove them here as well. 
        # Here we just need to update the features in each image to make it consistent
        all_feat <- Features(image_new)
        all_feat_to_remove <- all_feat[!(all_feat %in% c(intersecting_features, control_features, snvs_to_preserve))]
        features_to_keep <- all_feat[!(all_feat %in% all_feat_to_remove)]
        image_new <- subset(image_new, features = features_to_keep)
        sample_dot <- gsub("\\_",".",gsub("\\-",".",sampleID))
        image_name_sample <- paste0(image_name,".",sample_dot)
        obj_new[[image_name]] <- NULL
        obj_new[[image_name_sample]] <- image_new
    }
    DefaultAssay(obj_new) <- "Xenium"
    obj_new <- NormalizeData(obj_new, assay = "Xenium")
    obj_new <- FindVariableFeatures(obj_new, assay = "Xenium")
    obj_new <- ScaleData(obj_new, assay = "Xenium")
    obj_new <- SCTransform(obj_new, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2")
    obj_new <- FindVariableFeatures(obj_new, selection.method = "vst", nfeatures = dim(obj_new)[1], assay = "SCT")
    obj_new <- RunPCA(obj_new, npcs = opt$pc_num_individual, features = rownames(obj_new))
    obj_new <- RunUMAP(obj_new, dims = 1:opt$pc_num_individual, reduction.name = paste("umap.",opt$pc_num_individual,"PC", sep = ""), reduction.key = paste("UMAP",opt$pc_num_individual,"PC_",sep=""))
    obj_new <- FindNeighbors(obj_new, reduction = "pca", dims = 1:opt$pc_num_individual)
    obj_new <- FindClusters(obj_new, resolution = 0.3)
    DefaultAssay(obj_new) <- "SCT"
    Idents(obj_new) <- "seurat_clusters"
    print(paste0("validObject:"))
    print(validObject(obj_new))
    # This part can be un-commented if you want to check that the object is structured and plots correctly. You may need to change some of the genes/snvs so that at least one of them in each call is present in the seurat object.
    # pdf(paste(out_dir,"/",prefix,"_individual_processed_",sampleID,".pdf", sep=""), useDingbats=FALSE, height=8, width = 15)
    # print(rasterize(DimPlot(obj_new, reduction = paste("umap.",opt$pc_num_individual,"PC", sep=""), group.by = "seurat_clusters", label=TRUE, label.size=6, raster=FALSE), layers="Point", dpi=300))
    # print(rasterize(DimPlot(obj_new, reduction = paste("umap.",opt$pc_num_individual,"PC", sep=""), group.by = "sample_ID", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    # print(rasterize(FeaturePlot(obj_new, reduction = paste("umap.",opt$pc_num_individual,"PC", sep=""), features = c("KRT7","SOX9"), raster = F), layers="Point", dpi=300))
    # print(rasterize(FeaturePlot(obj_new, reduction = paste("umap.",opt$pc_num_individual,"PC", sep=""), features = c("MKI67"), raster = F), layers="Point", dpi=300))
    # assay_list <- Assays(obj_new) 
    # snv_assay <- assay_list[grepl("snv",assay_list)][1]
    # DefaultAssay(obj_new) <- snv_assay # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
    # image_list <- Images(obj_new)
    # snv_image <- image_list[grepl("snv",image_list)][1]
    # non_snv_image <- image_list[!(grepl("snv",image_list))][1]
    # print("testing plotting of non-SNV related information")
    # DefaultFOV(obj_new, assay="SCT") <- non_snv_image
    # DefaultBoundary(obj_new[[non_snv_image]]) <- "segmentation"
    # print(rasterize(ImageDimPlot(obj_new, fov = non_snv_image, group.by="cell_type_individual", border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600))
    # print(rasterize(ImageFeaturePlot(obj_new, fov = non_snv_image, features = c("KRAS","SOX9"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
    # print(rasterize(ImageFeaturePlot(obj_new, fov = non_snv_image, features = c("MKI67"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
    # print("testing plotting of non-SNV related information")
    # DefaultFOV(obj_new, assay=snv_assay) <- snv_image
    # DefaultBoundary(obj_new[[snv_image]]) <- "segmentation"
    # print(rasterize(ImageDimPlot(obj_new, fov = snv_image, group.by="cell_type_individual", border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600))
    # print(rasterize(ImageFeaturePlot(obj_new, fov = snv_image, features = c("KRAS","KRAS-p-G12V-WT","KRAS-p-G12D-ALT-21-A"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
    # new_colors = c("magenta","cyan","yellow")
    # new_markers = c("KRAS","KRAS-p-G12V-WT","KRAS-p-G12D-ALT-21-A")
    # print(ImageDimPlot(obj_new, fov = snv_image, #axes = TRUE, 
    #                 border.color = "white", border.size = 0.001, cols = rep("black", times = length(unique(obj_new$seurat_clusters))),
    #                 mols.size = 0.1,
    #                 molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
    #                 axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1))
    # dev.off()
    seurat.list[[i]] <- obj_new
    total_cells = total_cells + dim(obj_new)[2]
}

print(total_cells)
print("complete seurat.list")
print(seurat.list)
saveRDS(seurat.list, paste0(out_dir, "/",prefix,"_seurat.list_checkpoint.rds"))
print("Merging to one object")
# found out after initial test that JoinLayers only works on v5 Assays and the the individual objects are unintentionally a mix of v5 and v3.
# Error in UseMethod(generic = "JoinLayers", object = object) :
#   no applicable method for 'JoinLayers' applied to an object of class "c('Assay', 'KeyMixin')"
# It turns out that the Xenium assay is a Seurat3 assay instead of a Seurat5 assay so if 
# we merge/integrate the object is already combined into a single assay because v3 Seurat assays are not layers that can be kept separate
for (i in iterables) {
    obj_new <- seurat.list[[i]]
    assay_list <- Assays(obj_new)
    for (assay in assay_list) {
        assay_class = class(obj_new[[assay]])
        if ((assay_class == "Assay")) {
            obj_new[[assay]] <- as(object = obj_new[[assay]], Class = "Assay5")
        }
    }
    seurat.list[[i]] <- obj_new
}

print("Merging to one object")
obj_1 <- seurat.list[[1]]
obj.merge <- merge(obj_1, y = seurat.list[2:length(seurat.list)], add.cell.ids = sample_vector, merge.data = T, project = "PDAC_primary_snRNA")
rm(obj_1,seurat.list)
gc()
obj.merge$barcode <- rownames(obj.merge@meta.data)
print(SCTResults(obj.merge, slot="umi.assay"))
# $model1
# [1] "Xenium"
# 
# $model1.1
# [1] "Xenium"
# 
# $model1.2
# [1] "Xenium"
# 
# $model1.3
# [1] "Xenium"
# 
# $model1.4
# [1] "Xenium"
# 
# $model1.5
# [1] "Xenium"
# 
# $model1.6
# [1] "Xenium"
# 
# $model1.7
# [1] "Xenium"
# 
# $model1.8
# [1] "Xenium"
# after merging there is one model per dataset used in the merge.
# Joining the layers so that everything is treated like a single thing.
assay_list <- Assays(obj.merge)
for (assay in assay_list) {
    assay_class = class(obj.merge[[assay]])
    if (assay_class == "Assay5") {
        obj.merge[[assay]] <- JoinLayers(obj.merge[[assay]])
    }
}
print("SCTransform")
# lots of cells so setting conserve.memory = T
obj.merge <- SCTransform(obj.merge, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2", conserve.memory = T)
print(SCTResults(obj.merge, slot="umi.assay"))
# [1] "Xenium"
# Becuase we ran join layers there is only one SCT layer/matrix/model in the assay that will be used.
# obj.merge <- RunPCA(obj.merge, assay = "SCT", npcs = 200, verbose = TRUE)
# normally I check PC contribution up to npcs = 200, however, when I tried that here I got the following warning:
# Warning in irlba(A = t(x = object), nv = npcs, ...) :
# You're computing too large a percentage of total singular values, use a standard svd instead.
# In order to avoid this and keep it everything consistent between merged and non-merged analysis I'm using npcs = 100
# which does not result in the warning.
# dont need the following commented lines of code because there is a single SCT model in this object instead of many
# based on the recommendations in one of the comments here: https://github.com/satijalab/seurat/issues/8235 need to update the following if we want 
# for (i in iterables) {
#     slot(object = obj.merge@assays$SCT@SCTModel.list[[i]], name="umi.assay") <- "Xenium"
# }
# print(SCTResults(obj.merge, slot="umi.assay"))
# $`1`
# [1] "Xenium"
# $`2`
# [1] "Xenium"
# $`3`
# [1] "Xenium"
# $`4`
# [1] "Xenium"
# $`5`
# [1] "Xenium"
# $`6`
# [1] "Xenium"
# $`7`
# [1] "Xenium"
obj.merge <- RunPCA(obj.merge, assay = "SCT", npcs = 100, verbose = TRUE)
# this generates a single reduction though.
print("plotting")
pdf(paste0(out_dir,"/","Elbowplot_",prefix,"_single_SCT.pdf"))
print(ElbowPlot(obj.merge, ndims = 100, reduction = "pca"))
dev.off()
# obj.merge <- PrepSCTFindMarkers(obj.merge, assay = "SCT", verbose = TRUE)
# saveRDS(seurat.merge,paste0(out_dir, "/",prefix,"_no_umap.rds"))
clusters = c("5","10","15","20", "30", "40", "50", "60", "70", "80", "90", "100")
dims = c(5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)
for (i in 1:12) {
    dimen = dims[i]
    cluster = clusters[i]
    obj.merge <- RunUMAP(obj.merge, reduction = "pca", dims = 1:dimen, reduction.name = paste("umap.",cluster,"PC", sep = ""), reduction.key = paste("UMAP",cluster,"PC_",sep=""))
    obj.merge <- FindNeighbors(obj.merge, reduction = "pca", dims = 1:dimen, force.recalc = T)
    obj.merge <- FindClusters(obj.merge, resolution = 0.3)
    pdf(paste0(out_dir,"/Dimplots_",prefix,"_single_SCT_",cluster,"_PC.pdf"), useDingbats=FALSE, height=8, width = 15)
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "seurat_clusters", label=TRUE, label.size=6, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "sample_ID", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "piece_ID", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "case_ID", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "cell_type_individual", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "neoplasm_normal_unknown", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "neoplasm_normal_unknown_subclone", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "Site", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "Treatment", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "Progression", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "Sex", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "Age", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(DimPlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), group.by = "Section_KRAS_call", label=TRUE, label.size=2, raster=FALSE), layers="Point", dpi=300))
    print(rasterize(FeaturePlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nFeature_Xenium"), raster = F), layers="Point", dpi=300))
    print(rasterize(FeaturePlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nCount_Xenium"), raster = F), layers="Point", dpi=300))
    print(rasterize(FeaturePlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nFeature_ControlProbe"), raster = F, order = T), layers="Point", dpi=300))
    print(rasterize(FeaturePlot(obj.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nCount_ControlProbe"), raster = F, order = T), layers="Point", dpi=300))
    dev.off()
}
# not plotting these for now to save time.
# print(rasterize(FeaturePlot(seurat.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nFeature_ControlCodeword"), raster = F, order = T), layers="Point", dpi=300))
# print(rasterize(FeaturePlot(seurat.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nCount_ControlCodeword"), raster = F, order = T), layers="Point", dpi=300))
# print(rasterize(FeaturePlot(seurat.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nFeature_BlankCodeword"), raster = F, order = T), layers="Point", dpi=300))
# print(rasterize(FeaturePlot(seurat.merge, reduction = paste("umap.",cluster,"PC", sep=""), features = c("nCount_BlankCodeword"), raster = F, order = T), layers="Point", dpi=300))
obj.merge <- FindNeighbors(obj.merge, reduction = "pca", dims = 1:50)#, force.recalc = T) #seurat no longer supports
obj.merge <- FindClusters(obj.merge, resolution = 0.3)
saveRDS(obj.merge, paste0(out_dir,"/",prefix,"_single_SCT.rds")) #add umaps to the end
write.table(obj.merge@meta.data, paste0(out_dir,"/",prefix,"_single_SCT_20240205_metadata.tsv"), quote=F, sep='\t')
library(future)
options(future.globals.maxSize= +Inf)
# plan(multicore, workers = 40)
plan(sequential)
Idents(obj.merge) <- "seurat_clusters"
degs <- FindAllMarkers(obj.merge, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
write.table(degs,paste0(out_dir,"/",prefix,"_single_SCT_FindAllMarkers_seurat_clusters_50PC_res0.3_20241209.tsv"),sep="\t",quote=FALSE)
plan(sequential)
#This part can be un-commented if you want to check that the object is structured and plots images correctly post merging. You may need to change some of the genes/snvs so that at least one of them in each call is present in the seurat object.
pdf(paste(out_dir,"/",prefix,"_single_SCT_plotting_check.pdf", sep=""), useDingbats=FALSE, height=8, width = 15)
assay_list <- Assays(obj.merge)
snv_assay <- assay_list[grepl("snv",assay_list)][1]
DefaultAssay(obj.merge) <- snv_assay # this needs to be done before you crop or the crop will not be associated with the correct assay and you can't plot variant probes
image_list <- Images(obj.merge)
snv_image <- image_list[grepl("snv",image_list)] # there is one fov per initial object
non_snv_image <- image_list[!(grepl("snv",image_list))] #there is one fov per initial object
print("testing plotting of non-SNV related information")
for (image_name in non_snv_image) {
    DefaultFOV(obj.merge, assay="SCT") <- image_name
    DefaultBoundary(obj.merge[[image_name]]) <- "segmentation"
} 
print(rasterize(ImageDimPlot(obj.merge, fov = non_snv_image, group.by="cell_type_individual", border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600))
print(rasterize(ImageDimPlot(obj.merge, fov = non_snv_image, group.by="seurat_clusters", border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600))
# the next line produces a warning in my case because SOX9 is not present in the merged object. This warning can be ignored. If none of the genes passed to features = c() are present in the object then it will error.
# it also produces another warning because the snv_assay has not been normalized.
print(rasterize(ImageFeaturePlot(obj.merge, fov = non_snv_image, features = c("KRAS","SOX9"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
print(rasterize(ImageFeaturePlot(obj.merge, fov = non_snv_image, features = c("MKI67"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
DefaultAssay(obj.merge) <- "SCT"
print(rasterize(ImageFeaturePlot(obj.merge, fov = non_snv_image, features = c("KRAS","SOX9"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
print(rasterize(ImageFeaturePlot(obj.merge, fov = non_snv_image, features = c("MKI67"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
dev.off()
print("testing plotting of non-SNV related information")
pdf(paste(out_dir,"/",prefix,"_single_SCT_plotting_check_snvs.pdf", sep=""), useDingbats=FALSE, height=8, width = 15)
for (image_name in snv_image) {
    DefaultFOV(obj.merge, assay=snv_assay) <- image_name
    DefaultBoundary(obj.merge[[image_name]]) <- "segmentation"
} 
print(rasterize(ImageDimPlot(obj.merge, fov = snv_image, group.by="cell_type_individual", border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600))
print(rasterize(ImageFeaturePlot(obj.merge, fov = snv_image, features = c("KRAS","KRAS-p-G12V-WT","KRAS-p-G12D-ALT-21-A","KRAS-p-G12D-ALT-T"), min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
new_colors = c("magenta","cyan","yellow")
new_markers = c("KRAS","KRAS-p-G12V-WT","KRAS-p-G12D-ALT-21-A","KRAS-p-G12D-ALT-T")
# molecules appear to only be plotted according to the coordinates of the first FOV that is plotted. So if you try to plot all of the FOVs at once then the images for the later FOVs will show the molecules from the first FOV but not their own.
# If you plot the FOVs separately in something like a for loop then there is no issue.
# this follows the recommended solution from https://github.com/satijalab/seurat/issues/8713 
print(ImageDimPlot(obj.merge, fov = snv_image, #axes = TRUE,
                   border.color = "white", border.size = 0.001, cols = rep("black", times = length(unique(obj.merge$seurat_clusters))),
                   mols.size = 0.1,
                   molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
                   axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1))
for (image_name in snv_image) {
    print(ImageDimPlot(obj.merge, fov = image_name, #axes = TRUE,
                       border.color = "white", border.size = 0.001, cols = rep("black", times = length(unique(obj.merge$seurat_clusters))),
                       mols.size = 0.1,
                       molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
                       axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1))
}
dev.off()
