#!/usr/bin/env Rscript
library(tidyverse)
library(Seurat)
#library(Signac)
library(tools)
library(ggplot2)
set.seed(1234)
args <- commandArgs(trailingOnly = TRUE)
#cell_typing_path <- "/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_v12062021.txt"
#cell_typing_path <- "/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_caf_v20230508.txt"
#cell_typing_path <- "/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt"
#cell_typing_path <- "/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_Liver_v20220410.txt"
#cell_typing_path <- "/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_Lung_v20220412.txt"
#cell_typing_path <- "/diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_lineage_v20220523.txt"
cell_typing_path = args[1]
cell_typing_table <- read.table(cell_typing_path, sep = "\t", header = T, row.names=NULL)
cell_types_vector <- unique(cell_typing_table$Gene_set)
print(head(cell_typing_table))
#dir_with_all_objects <- args[1]
out_dir = args[2]
obj_path = args[3]
annotations_to_use_path <- args[4] #this needs to be a tsv file with a header
obj_name <- file_path_sans_ext(basename(obj_path))
cell_state_name <- file_path_sans_ext(basename(cell_typing_path))
prefix <- args[5]
print(obj_name)
if (file_ext(obj_path) == "rds") {
    obj <- readRDS(obj_path)
    if (!(is.null(obj@assays$SCT))) {
        plots <- NA
        DefaultAssay(obj) <- "SCT"
        print(annotations_to_use_path)
        if (file_ext(annotations_to_use_path) == "tsv") {
            annotations = read.table(annotations_to_use_path, sep = '\t', header = T)
            colnames(annotations) <- c('barcode', 'marker_ident')
            rownames(annotations) <- annotations$barcode
            annotations$barcode <- NULL
            obj <- AddMetaData(obj, annotations)
        } else if (file_ext(annotations_to_use_path) == "csv") {
            annotations = read.table(annotations_to_use_path, sep = ',', header = T)
            colnames(annotations) <- c('barcode', 'marker_ident')
            rownames(annotations) <- annotations$barcode
            annotations$barcode <- NULL
            obj <- AddMetaData(obj, annotations)
        }
        Idents(obj) <- "marker_ident"
        #Idents(obj) <- "sub.4.12.cluster"
        #Idents(obj) <- "sub.6.7.cluster"
        #Idents(obj) <- "sub.2.cluster"
        #Idents(obj) <- "sub.1.6.cluster"
        #Idents(obj) <- "sub.4.6.8.12.14.cluster"
        #Idents(obj) <- "sub.14.cluster"
        #Idents(obj) <- "sub.5.cluster"
        #Idents(obj) <- "sub.4.cluster"
        #Idents(obj) <- "sub.1.9.13.cluster"
        sample_genes <- rownames(x = obj)
        pdf(file = paste0(out_dir,"/",prefix,".",obj_name,".",cell_state_name,".pdf", sep = ""), useDingbats=FALSE,width=7,height=11)
        for (cell_type in cell_types_vector) {
            print(cell_type)
            genes_to_plot <- unique(cell_typing_table$Gene[cell_typing_table$Gene_set == cell_type])
            count <- 0
            for (gene in genes_to_plot) {
                if (gene %in% sample_genes) {
                    count = count + 1
                }
            }
            if (count > 0) {
                print(genes_to_plot)
                plot <- DotPlot(object = obj, assay='SCT', features = genes_to_plot, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100) + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red") + scale_x_discrete(guide = guide_axis(n.dodge=2)) + labs(title=paste(cell_type," markers in ",prefix,sep="")) + RotatedAxis()
                #plot <- DotPlot(object = obj, assay='SCT', features = genes_to_plot) + labs(title=paste(cell_type," markers in ",sample,sep="")) + RotatedAxis()
                print(plot)
            }
        }
        dev.off()
    }
}

# When you see this error it means that there is a single gene that is listed twice for the same cell type. 
# You need to go into the reference table that you are using and delete whatever line is the duplicate.
#[1] "Effector memory"
#Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final data frame: CX3CR1, SELL, IL2
#Error in `levels<-`(`*tmp*`, value = as.character(levels)) : 
#    factor level [9] is duplicated
#Calls: DotPlot -> factor
#In addition: There were 26 warnings (use warnings() to see them)
#Execution halted


# When you see this error it means that there are NA values in the vector that is currently assigned to the Idents(obj) part of the Seurat object
# Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final data frame: KRT19, KRT8, KRT18, KRT17, KRT7, KRT5, KRT6A, KRT14, EPCAM, TACSTD2, ANXA2, S100A10, S100A11, S100A16, TPM1, TFF1, S100A6, AGR2, C19orf33                        
# Error in FetchData(object = object, vars = features, cells = cells) :                              
#   None of the requested variables were found (10 out of 19 shown): KRT19, KRT8, KRT18, KRT17, KRT7, KRT5, KRT6A, KRT14, EPCAM, TACSTD2
# Calls: DotPlot -> FetchData                                                                        
# In addition: Warning message:                                                                      
# In cells.idents["NA"] <- names(x = which(x = is.na(x = Idents(object = object)[cells]))) :         
#   number of items to replace is not a multiple of replacement length                               
# Execution halted 
# To fix this you need to reset the Idents(obj) to something else or replace the NA values with other things.
# Example: Idents(obj) <- "seurat_clusters"
