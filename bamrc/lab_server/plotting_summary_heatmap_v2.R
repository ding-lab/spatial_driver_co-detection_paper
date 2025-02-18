library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(grid)
set.seed(1234)
vaf_table_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/bamrc_vaf_summary_table.tsv"
vaf_df = as.data.frame(read.table(vaf_table_path, sep="\t", header=F))
vaf_df = as.data.frame(t(vaf_df))
colnames(vaf_df) <- vaf_df[1,]
var_names <- unlist(vaf_df$index[2:length(vaf_df[,1])])
rownames(vaf_df) <- vaf_df$index
vaf_df$index <- NULL
vaf_df <- vaf_df[2:length(vaf_df$`HT227P1-S1H1L1U1`),]
for (col_name in colnames(vaf_df)) {
    vaf_df[,col_name][vaf_df[,col_name] == ""] <- NA
    vaf_df[,col_name] <- as.numeric(vaf_df[,col_name])
}
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/test/"
var_anno_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/variant_annotation.tsv"
sample_anno_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/cohort_annotation.tsv"
sample_summary_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/cohort_summary.tsv"
sample_summary_df <- read.table(sample_summary_path, sep ='\t', header = T, row.names = "Sample_ID")
sample_summary_df <- as.data.frame(t(sample_summary_df))
sample_summary_colors <- colorRamp2(c(FALSE,TRUE), c("grey80", "grey0"))
vaf_colors = colorRamp2(c(0,0.4,0.75,1), c("midnightblue","magenta","#FC8D62","gold"))
sample_anno_df <- as.data.frame(read.table(sample_anno_path, sep = "\t",header = T,row.names = "Sample_ID"))
set3_colors = brewer.pal(n=10,name = "Set3")
c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
set1_colors = brewer.pal(n=9,name = "Set1")
c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
set2_colors = brewer.pal(n=8,name = "Set2")
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
var_anno_df <- read.table(var_anno_path, sep ='\t', header = T, row.names = "Alternate.probe.name")
sample_order <- colnames(sample_summary_df)
variant_order <- rownames(vaf_df)
Spectral_colors = brewer.pal(n=11,name = "Spectral")
c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598","#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
sample_anno_df = sample_anno_df[sample_order,] # set the order of the annotations so the vector ordering is correct on the heatmap
sample_column_anno = HeatmapAnnotation(Site = sample_anno_df$Site,
                                       Type = sample_anno_df$Type,
                                       Disease = sample_anno_df$Disease,
                                       Sample_control = sample_anno_df$Sample_control,
                                       col = list(Site = c("Breast" = "#BC80BD",
                                                           "Liver" = "#8DD3C7",
                                                           "Kidney" = "#FFFFB3",
                                                           "Colon" = "#FDB462",
                                                           "Brain" = "#D9D9D9",
                                                           "Bone Marrow" = "#FB8072",
                                                           "Scalp Plasmocytoma" = "#BEBADA",
                                                           "Pancreas" = "#B3DE69",
                                                           "Body wall" = "#80B1D3",
                                                           "Bone Marrow hip" = "#FCCDE5"),
                                                  Type = c("Primary" = "#0075BF",
                                                           "Metastasis" = "#ED1652",
                                                           "Normal" = "#A3D39C"),
                                                  Disease = c("BRCA" = "#F781BF",
                                                              "ccRCC" = "#FFFF33",
                                                              "CHOL" = "#999999",
                                                              "CRC" = "#A65628",
                                                              "GBM" = "#984EA3",
                                                              "MM" = "#E41A1C",
                                                              "PDAC" = "#377EB8",
                                                              "Healthy" = "#4DAF4A"),
                                                  Sample_control = c("Cancer sample" = '#CC79A7',
                                                                     "Healthy control" = '#009E73')),
                                       gp = gpar(col = "grey30"),
                                       na_col = "black",
                                       annotation_name_side = "left")
var_anno_df = var_anno_df[variant_order,] # set the order of the annotations so the vector ordering is correct on the heatmap
var_row_anno = rowAnnotation(Variant_type_name = var_anno_df$Variant_type_name,
                             Mutation_role = var_anno_df$Mutation_role,
                             col = list(Variant_type_name = c("SNV" = "#E6F598",
                                                              "Indel" = "#5E4FA2"),
                                        Mutation_role = c("Driver" = "#FD63CB",
                                                          "Germline" = "#22FDEB",
                                                          "Clonal" = "#FEE100")),
                             gp = gpar(col = "grey30"),
                             na_col = "black")
colorspace::diverge_hcl(12,h=c(180,70),c=70,l=c(90,95))
c("#22FDEB", "#83FAEC", "#AEF7EE", "#CBF5EF", "#E0F3F0", "#EDF1F0", "#F1F0EE", "#F3EEE3", "#F6ECD2", "#F9E9BD", "#FBE5A2", "#FEE100")
colorspace::diverge_hcl(12,h=c(128,330),c=98,l=c(65,90))
c("#10B717", "#64C366", "#90CE91", "#B2D6B2", "#CCDDCC", "#DEE1DE", "#E4DFE2", "#E9D2DF", "#EFC0DB", "#F4A8D6", "#F98AD1", "#FD63CB")
# sample_summary_colors <- structure(0:1, names = c("False", "True")) # black, red, green, blue
vaf_df = vaf_df[,sample_order]
sample_summmary = Heatmap(as.matrix(sample_summary_df), 
                          name = "Data Availability", 
                          col = sample_summary_colors,
                          show_column_dend = FALSE,
                          show_row_dend = FALSE,
                          row_order = c("Xenium_snvs","bulk_WES","snRNAseq","scRNAseq"),
                          column_order = sample_order,
                          row_names_side = "left",
                          heatmap_legend_param = list(title = "Data availability", color_bar = "discrete", at = c(0,1)),
                          rect_gp = gpar(col = "grey30", lwd = 1),
                          column_names_gp = gpar(fontsize = 8),
                          row_names_gp = gpar(fontsize = 8))
vaf_summary = Heatmap(as.matrix(vaf_df), 
                      name = "Bulk WES VAF", 
                      col = vaf_colors,
                      na_col = "grey80",
                      show_column_dend = FALSE,
                      show_row_dend = FALSE,
                      row_order = variant_order,
                      column_order = sample_order,
                      row_names_side = "left",
                      rect_gp = gpar(col = "grey30", lwd = 1),
                      column_names_gp = gpar(fontsize = 8),
                      row_names_gp = gpar(fontsize = 8))
sample_summmary_2 = Heatmap(as.matrix(sample_summary_df), 
                            name = "Data Availability", 
                            col = sample_summary_colors,
                            show_column_dend = FALSE,
                            show_row_dend = FALSE,
                            row_order = c("Xenium_snvs","bulk_WES","snRNAseq","scRNAseq"),
                            column_order = sample_order,
                            row_names_side = "left",
                            top_annotation = sample_column_anno,
                            heatmap_legend_param = list(title = "Data availability", color_bar = "discrete", at = c(0,1)),
                            rect_gp = gpar(col = "grey30", lwd = 1),
                            column_names_gp = gpar(fontsize = 8),
                            row_names_gp = gpar(fontsize = 8))
vaf_summary_2 = Heatmap(as.matrix(vaf_df), 
                        name = "Bulk WES VAF", 
                        col = vaf_colors,
                        na_col = "grey80",
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        row_order = variant_order,
                        column_order = sample_order,
                        row_names_side = "left",
                        left_annotation = var_row_anno,
                        rect_gp = gpar(col = "grey30", lwd = 1),
                        column_names_gp = gpar(fontsize = 8),
                        row_names_gp = gpar(fontsize = 8))
pdf(paste0(out_dir,"/","Cohort_summary_heatmap_v3.pdf"),height = 10, width = 13)
print(sample_summmary)
print(vaf_summary)
print(sample_summmary_2)
print(vaf_summary_2)
print(sample_summmary_2 %v% vaf_summary_2)
dev.off()
# library(circlize)
library(grid)
vaf_min <- min(vaf_df[vaf_df>0 & !(is.na(vaf_df))])
print(vaf_min)
effective_zero = vaf_min/100
print(effective_zero)
# highlight_matrix <- vaf_df > 0
# highlight_matrix[is.na(vaf_df)] <- FALSE
# highlight_matrix <- as.matrix(highlight_matrix)
offset = effective_zero*1e-10
print(offset)
vaf_colors = colorRamp2(c(0, effective_zero, effective_zero+offset, 0.4,0.75, 1), c("white","white","midnightblue","magenta","#F69795","gold"))
vaf_summary_2 = Heatmap(as.matrix(vaf_df), 
                        name = "Bulk WES VAF", 
                        col = vaf_colors,
                        # for some reason the below is not working and I do not think I have time to troubleshoot it right now.
                        # it is from this post: https://stackoverflow.com/questions/78696356/complexheatmap-highlight-cells-by-changing-alpha 
                        # cell_fun = function(j, i, x, y, width, height, fill) {
                        #     if (!highlight_matrix[i, j]) {
                        #         fill <- circlize::add_transparency(fill, 0.1)
                        #     }
                        #     grid::grid.rect(x, y,
                        #                     width = width, height = height,
                        #                     gp = gpar(fill = fill)
                        #     )
                        # }
                        na_col = "grey80",
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        row_order = variant_order,
                        column_order = sample_order,
                        row_names_side = "left",
                        left_annotation = var_row_anno,
                        rect_gp = gpar(col = "grey30", lwd = 1))
pdf(paste0(out_dir,"/","Cohort_summary_heatmap_v4.pdf"),height = 10, width = 13)
print(sample_summmary)
print(vaf_summary)
print(sample_summmary_2)
print(vaf_summary_2)
print(sample_summmary_2 %v% vaf_summary_2) #this is the one used for figure 1b
dev.off()
