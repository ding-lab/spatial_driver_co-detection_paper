# figure 1c
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
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/"
var_anno_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/variant_annotation.tsv"
sample_anno_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/cohort_annotation.tsv"
sample_summary_path = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bamrc_output/vaf_summary_table/v3_freeze/cohort_summary.tsv"
sample_summary_df <- read.table(sample_summary_path, sep ='\t', header = T, row.names = "Sample_ID")
sample_summary_df <- as.data.frame(t(sample_summary_df))
sample_summary_colors <- colorRamp2(c(FALSE,TRUE), c("grey80", "grey0"))
vaf_colors = colorRamp2(c(0,0.5, 1), c("midnightblue","magenta","gold"))
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
                          rect_gp = gpar(col = "grey30", lwd = 1))
vaf_summary = Heatmap(as.matrix(vaf_df), 
                      name = "Bulk WES VAF", 
                      col = vaf_colors,
                      na_col = "grey80",
                      show_column_dend = FALSE,
                      show_row_dend = FALSE,
                      row_order = variant_order,
                      column_order = sample_order,
                      row_names_side = "left",
                      rect_gp = gpar(col = "grey30", lwd = 1))
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
                            rect_gp = gpar(col = "grey30", lwd = 1))
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
                        rect_gp = gpar(col = "grey30", lwd = 1))
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
vaf_colors = colorRamp2(c(0, effective_zero, effective_zero+offset, 0.5, 1), c("white","white","midnightblue","magenta","gold")) #"transparency" workaround because I am trying to make all of the 0 values transparent is based off of this post: https://stackoverflow.com/questions/59891243/complexheatmap-highlight-specific-values-in-a-heatmap
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
pdf(paste0(out_dir,"/","Cohort_summary_heatmap_v4.pdf"),height = 10, width = 13) #used this one.
print(sample_summmary)
print(vaf_summary)
print(sample_summmary_2)
print(vaf_summary_2)
print(sample_summmary_2 %v% vaf_summary_2)
dev.off()

# figure 1 supplement violin plots and bar plots of detection rate:
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
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/"
prefix = "All_samples_merge_snRNA_20250119_separateSCT"
obj_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/"
obj.merge = readRDS(paste0(obj_dir,"/",prefix,"_separate_SCT_umaps.rds"))
Idents(obj.merge) <- "neoplasm_normal_unknown"
obj.merge$neoplasm_normal_unknown <- factor(obj.merge$neoplasm_normal_unknown, levels = c("neoplastic","normal"))
library(scales)
cols <- c("neoplastic" = hue_pal()(2)[1], "normal" = hue_pal()(2)[2])
# neoplastic     normal 
# "#F8766D"  "#00BFC4" # red  blue
genes_to_plot_uni <- unique(c('APC','APC','ATM','ATM','BAP1','BRAF','BRCA2','CTNNB1','DIS3','DIS3','DLGAP4','EEF1A1','FANCA','GNAS','IDH1','IDH1','IRF4','KRAS','KRAS','KRAS','LDHB','MAP1B','NRAS','PELI1','PGAP2','PIK3CA','SMG1','TP53','TP53','TPD52L1'))
p1 <- VlnPlot(obj.merge, features = genes_to_plot_uni, assay = "SCT", cols = cols, group.by = "sample_ID", split.by = "neoplasm_normal_unknown", y.max=4, pt.size=0, raster=FALSE)
pdf(paste0(out_dir,"/",prefix,"_VlnPlot_variant_targeted_genes_by_case.pdf"),useDingbats = F, height=35,width=50)
print(p1)
dev.off()

# figure 1 supplement violin plots and bar plots of detection rate vs expression level
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/
# conda activate seurat5
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(Seurat)
library(optparse)
library(patchwork)
library(ggpmisc)
set.seed(1234)
library(future)
library(scales)
library(ggrepel)
library(ggh4x)
options(future.globals.maxSize= +Inf)
plan(sequential)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/"
all_sample_summary = read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header = T)
variant_xenium_snRNA = read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/sample_mutation_table_based_on_bulk_7_with_snRNA.tsv", sep = '\t', header=T)
snRNA_to_xenium_linking = read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/snRNA_to_xenium_linking_table.tsv", sep='\t', header=T)
obj_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/"
snRNA_prefix = "All_samples_merge_snRNA_20250119_separateSCT"
out_dir =  "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/absolute_expression/"
snRNA_obj = readRDS(paste0(obj_dir,"/",snRNA_prefix,"_separate_SCT_umaps.rds"))
dir.create(out_dir)
variants_with_snRNA = variant_xenium_snRNA$Variant[!(variant_xenium_snRNA$snRNA_sample == "")]
# the cell types for CPT0079400003 were apparently the old ones where the rhabdoid cancer cell clone was not appropriately labeled.
snRNA_obj$cell_type_individual[(snRNA_obj$seurat_clusters_individual == 8) & (snRNA_obj$cell_type_individual == "ccRCC cancer cell") & (snRNA_obj$sample_ID == "CPT0079400003")] <- "rhabdoid cancer cell"
snRNA_obj$neoplasm_normal_unknown[(snRNA_obj$seurat_clusters_individual == 8) & (snRNA_obj$cell_type_individual == "rhabdoid cancer cell") & (snRNA_obj$sample_ID == "CPT0079400003")] <- "neoplastic"
snRNA_obj$neoplasm_normal_unknown <- factor(snRNA_obj$neoplasm_normal_unknown, levels = c("neoplastic","normal"))
cols <- c("neoplastic" = hue_pal()(2)[1], "normal" = hue_pal()(2)[2])
# "#ED1652" "#F8766D"  "#00BFC4" "#0075BF"
Idents(snRNA_obj) <- "sample_ID"
all_sample_summary$absolute_proportion.cancer.cells.with.the.variant.probe <- all_sample_summary$Number.of.cancer.cells.with.the.variant.probe/all_sample_summary$total.number.of.cancer.cells.in.sample
all_sample_summary$absolute_proportion.cancer.cells.with.the.reference.probe <- all_sample_summary$Number.of.cancer.cells.with.the.reference.probe/all_sample_summary$total.number.of.cancer.cells.in.sample
all_sample_summary$absolute_proportion.non.cancer.cells.with.the.variant.probe <- all_sample_summary$Number.of.non.cancer.cells.with.the.variant.probe/all_sample_summary$total.number.of.normal.cells.in.sample
all_sample_summary$absolute_proportion.non.cancer.cells.with.the.reference.probe <- all_sample_summary$Number.of.non.cancer.cells.with.the.reference.probe/all_sample_summary$total.number.of.normal.cells.in.sample
pdf(paste0(out_dir,"/","Expression_vs_detection_",snRNA_prefix,"_VlnPlot_bar.pdf"),useDingbats = F,height=10, width=10)
for(variant in variants_with_snRNA) {
    print(variant)
    snRNA_samples = unlist(strsplit(variant_xenium_snRNA$snRNA_sample[variant_xenium_snRNA$Variant == variant], ","))
    plot_xen_samples = snRNA_to_xenium_linking$Sample_ID[snRNA_to_xenium_linking$matching_snRNA_sample_ID_1_closest_match %in% snRNA_samples]
    variant_sample_summary <- all_sample_summary[(all_sample_summary$Variant.probe.name == variant) & (all_sample_summary$sample.ID %in% plot_xen_samples),]
    number_xenium_plots = length(plot_xen_samples)
    number_of_samples <- length(unique(variant_sample_summary$sample.ID))
    variant_plotting_df <- data.frame(absolute_proportion = c(variant_sample_summary$absolute_proportion.cancer.cells.with.the.variant.probe, variant_sample_summary$absolute_proportion.cancer.cells.with.the.reference.probe, variant_sample_summary$absolute_proportion.non.cancer.cells.with.the.variant.probe, variant_sample_summary$absolute_proportion.non.cancer.cells.with.the.reference.probe),
                                      sample.ID = c(variant_sample_summary$sample.ID, variant_sample_summary$sample.ID, variant_sample_summary$sample.ID, variant_sample_summary$sample.ID),
                                      group = rep(c("cancer with var probe","cancer with ref probe","normal with var probe","normal with ref probe"), each = number_of_samples))
    variant_plotting_df$group <- factor(variant_plotting_df$group, levels = c("cancer with var probe","cancer with ref probe","normal with var probe","normal with ref probe"))
    xenium_plotting_order <- c()
    for (snRNA_sample in snRNA_samples) {
        xenium_plotting_order <- c(xenium_plotting_order, snRNA_to_xenium_linking$Sample_ID[snRNA_to_xenium_linking$matching_snRNA_sample_ID_1_closest_match == snRNA_sample])
    }
    print(xenium_plotting_order)
    variant_plotting_df$sample.ID <- factor(variant_plotting_df$sample.ID, levels = xenium_plotting_order)
    gene <- unlist(strsplit(variant, "-"))[1]
    print(snRNA_samples)
    p1 <- VlnPlot(snRNA_obj, features = gene, assay = "SCT", cols = cols, idents = snRNA_samples, split.by = "neoplasm_normal_unknown", y.max=4, pt.size=0, raster=FALSE)
    p2 <- ggplot(variant_plotting_df, aes(y = absolute_proportion, x = sample.ID, fill = group)) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_manual(values = c("cancer with var probe" = "#ED1652",
                              "cancer with ref probe"= "#F8766D",
                              "normal with var probe" = "#0075BF",
                              "normal with ref probe" = "#00BFC4")) +
        theme_classic() +
        theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(color = "black", size = 12),
              #legend.position="none",
              axis.title=element_text(size=12),
              axis.ticks = element_line(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
        labs( y = paste("Propotion of all cells",sep = ""),
              #x = "Treatment", 
              title = paste0("Rate of ",variant," detection in cancer and normal cells"))
    list_of_plots <- list(p1)
    list_of_plots[[ length(list_of_plots) + 1 ]] <- p2
    # pdf(paste0(out_dir,"/","test",".pdf"),useDingbats = F,height=10, width=10)
    # print(p1)
    # print(p2)
    print(ggarrange(plotlist = list_of_plots, ncol = 1, nrow = length(list_of_plots)))
    # dev.off()
}
dev.off()
neoplastic_snRNA <- subset(snRNA_obj, subset = neoplasm_normal_unknown == "neoplastic")
normal_snRNA <- subset(snRNA_obj, subset = neoplasm_normal_unknown == "normal")
group_ident = "sample_ID"
all_sample_summary$proportion.of.cancer.with.variant.site.detected <- all_sample_summary$total.number.of.cancer.cells.with.either.variant.site.probe/all_sample_summary$total.number.of.cancer.cells.in.sample
all_sample_summary$proportion.of.normal.with.variant.site.detected <- all_sample_summary$total.number.of.normal.cells.with.either.variant.site.probe/all_sample_summary$total.number.of.normal.cells.in.sample
all_xen_samples_with_snRNA <- snRNA_to_xenium_linking$Sample_ID[!(snRNA_to_xenium_linking$matching_snRNA_sample_ID_1_closest_match == "")]
all_sample_summary$snRNA_sample <- NA
all_sample_summary$snRNA_sample_hyphens <- NA
all_sample_summary$sample.ID_hyphens <- NA
for (xen_sample in all_xen_samples_with_snRNA) {
    snRNA_sample_ID <- snRNA_to_xenium_linking$matching_snRNA_sample_ID_1_closest_match[snRNA_to_xenium_linking$Sample_ID == xen_sample]
    snRNA_sample_ID_hyphens <- gsub("_", "-", snRNA_sample_ID)
    xen_sample_hyphens <- gsub("_","-",xen_sample)
    all_sample_summary$snRNA_sample[all_sample_summary$sample.ID == xen_sample] <- snRNA_sample_ID
    all_sample_summary$snRNA_sample_hyphens[all_sample_summary$sample.ID == xen_sample] <- snRNA_sample_ID_hyphens
    all_sample_summary$sample.ID_hyphens[all_sample_summary$sample.ID == xen_sample] <- xen_sample_hyphens
}
sorted_all_xen_samples_with_snRNA <- sort(all_xen_samples_with_snRNA)
sorted_snRNA_samples <- sort(unique(all_sample_summary$snRNA_sample_hyphens)[!(is.na(unique(all_sample_summary$snRNA_sample_hyphens)))])
# variant = variants_with_snRNA[12]
pdf(paste0(out_dir,"/","Expression_log10_vs_total_detection_",snRNA_prefix,"_scatterplot.pdf"),useDingbats = F,height=10, width=10)
for (variant in variants_with_snRNA) {
    print(variant)
    gene <- unlist(strsplit(variant, "-"))[1]
    mutant_snRNA_samples = unlist(strsplit(variant_xenium_snRNA$snRNA_sample[variant_xenium_snRNA$Variant == variant], ","))
    mutant_snRNA_samples_hyphens <- gsub("_", "-", mutant_snRNA_samples)
    neoplastic_variant_sample_list <- AggregateExpression(neoplastic_snRNA, assays = "SCT", features = gene, group.by=group_ident)
    normal_variant_sample_list <- AggregateExpression(normal_snRNA, assays = "SCT", features = gene, group.by=group_ident)
    neoplastic_gene_sample <- neoplastic_variant_sample_list$SCT
    normal_gene_sample <- normal_variant_sample_list$SCT
    variant_sample_summary <- all_sample_summary[(all_sample_summary$Variant.probe.name == variant) & (all_sample_summary$sample.ID %in% all_xen_samples_with_snRNA),]
    rownames(neoplastic_gene_sample) <- gene
    rownames(normal_gene_sample) <- gene
    sort_neoplastic_gene_sample <- neoplastic_gene_sample[, sorted_snRNA_samples]
    sort_normal_gene_sample <- normal_gene_sample[gene, sorted_snRNA_samples]
    neoplastic_site_proportion <- variant_sample_summary[match(sorted_snRNA_samples, variant_sample_summary$snRNA_sample_hyphens) ,"proportion.of.cancer.with.variant.site.detected"]
    normal_site_proportion <- variant_sample_summary[match(sorted_snRNA_samples, variant_sample_summary$snRNA_sample_hyphens) ,"proportion.of.normal.with.variant.site.detected"]
    variant_site_plotting_df <- data.frame(gene_expression_log10 = log10(c(sort_neoplastic_gene_sample,sort_normal_gene_sample)),
                                      cell_label = rep(c("cancer cell", "normal cell"), each=length(sorted_snRNA_samples)),
                                      variant_site_detection_rate = c(neoplastic_site_proportion,normal_site_proportion),
                                      sample_ID = c(sorted_snRNA_samples,sorted_snRNA_samples))
    variant_site_plotting_df$cell_label = factor(variant_site_plotting_df$cell_label, levels = c("cancer cell", "normal cell")) 
    p1 <- ggplot(variant_site_plotting_df, aes(x = gene_expression_log10, y = variant_site_detection_rate, fill == cell_label, colour = cell_label, label = sample_ID)) +
        geom_point(size = 4) + 
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2","P")) , rsquared.conf.level = NA) +
        theme(axis.text.x = element_text(color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              axis.title=element_text(size=12),
              axis.ticks = element_line(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=2))+
        labs( y = paste("proportion of all cells in sample", sep = ""),
              x = "log10FC snRNA aggregate gene expression", 
              title = paste0(variant," Xenium site detection to snRNA gene expression correlation")) + 
        scale_fill_manual(values = c("cancer cell" = "#F8766D",
                                    "normal cell" = "#00BFC4")) +
        scale_colour_manual(values = c("cancer cell" = "#F8766D",
                                     "normal cell" = "#00BFC4")) + 
        geom_text_repel(data          = subset(variant_site_plotting_df, sample_ID %in% mutant_snRNA_samples_hyphens),
                        box.padding   = 1.5,
                        point.padding = 0.5,
                        force         = 150,
                        segment.size  = 0.2,
                        segment.color = "grey50") +
        force_panelsizes(rows = unit(7, "in"),
                         cols = unit(7, "in"))
    print(p1)
}
dev.off()
pdf(paste0(out_dir,"/","Expression_vs_total_detection_",snRNA_prefix,"_scatterplot.pdf"),useDingbats = F,height=10, width=10)
for (variant in variants_with_snRNA) {
    print(variant)
    gene <- unlist(strsplit(variant, "-"))[1]
    mutant_snRNA_samples = unlist(strsplit(variant_xenium_snRNA$snRNA_sample[variant_xenium_snRNA$Variant == variant], ","))
    mutant_snRNA_samples_hyphens <- gsub("_", "-", mutant_snRNA_samples)
    neoplastic_variant_sample_list <- AggregateExpression(neoplastic_snRNA, assays = "SCT", features = gene, group.by=group_ident)
    normal_variant_sample_list <- AggregateExpression(normal_snRNA, assays = "SCT", features = gene, group.by=group_ident)
    neoplastic_gene_sample <- neoplastic_variant_sample_list$SCT
    normal_gene_sample <- normal_variant_sample_list$SCT
    variant_sample_summary <- all_sample_summary[(all_sample_summary$Variant.probe.name == variant) & (all_sample_summary$sample.ID %in% all_xen_samples_with_snRNA),]
    rownames(neoplastic_gene_sample) <- gene
    rownames(normal_gene_sample) <- gene
    sort_neoplastic_gene_sample <- neoplastic_gene_sample[, sorted_snRNA_samples]
    sort_normal_gene_sample <- normal_gene_sample[gene, sorted_snRNA_samples]
    neoplastic_site_proportion <- variant_sample_summary[match(sorted_snRNA_samples, variant_sample_summary$snRNA_sample_hyphens) ,"proportion.of.cancer.with.variant.site.detected"]
    normal_site_proportion <- variant_sample_summary[match(sorted_snRNA_samples, variant_sample_summary$snRNA_sample_hyphens) ,"proportion.of.normal.with.variant.site.detected"]
    variant_site_plotting_df <- data.frame(gene_expression = c(sort_neoplastic_gene_sample,sort_normal_gene_sample),
                                           cell_label = rep(c("cancer cell", "normal cell"), each=length(sorted_snRNA_samples)),
                                           variant_site_detection_rate = c(neoplastic_site_proportion,normal_site_proportion),
                                           sample_ID = c(sorted_snRNA_samples,sorted_snRNA_samples))
    variant_site_plotting_df$cell_label = factor(variant_site_plotting_df$cell_label, levels = c("cancer cell", "normal cell")) 
    p1 <- ggplot(variant_site_plotting_df, aes(x = gene_expression, y = variant_site_detection_rate, fill == cell_label, colour = cell_label, label = sample_ID)) +
        geom_point(size = 4) + 
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2","P")) , rsquared.conf.level = NA) +
        theme(axis.text.x = element_text(color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              axis.title=element_text(size=12),
              axis.ticks = element_line(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=2))+
        labs( y = paste("proportion of all cells in sample", sep = ""),
              x = "snRNA aggregate gene expression", 
              title = paste0(variant," Xenium site detection to snRNA gene expression correlation")) + 
        scale_fill_manual(values = c("cancer cell" = "#F8766D",
                                     "normal cell" = "#00BFC4")) +
        scale_colour_manual(values = c("cancer cell" = "#F8766D",
                                       "normal cell" = "#00BFC4")) + 
        geom_text_repel(data          = subset(variant_site_plotting_df, sample_ID %in% mutant_snRNA_samples_hyphens),
                        box.padding   = 1.5,
                        point.padding = 0.5,
                        force         = 150,
                        segment.size  = 0.2,
                        segment.color = "grey50") +
        force_panelsizes(rows = unit(7, "in"),
                         cols = unit(7, "in"))
    print(p1)
}
dev.off()
# writing the all of the cell types the sn/scRNA samples to a single csv file
sample_vector <- unique(snRNA_obj$sample_ID)
snRNA_obj <- RenameCells(snRNA_obj, new.names =snRNA_obj$original_object_barcode)
snRNA_obj$barcode <- row.names(snRNA_obj@meta.data)
write.table(snRNA_obj@meta.data, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/All_samples_merge_snRNA_20250119_separateSCT_separate_SCT_20250125_metadata.tsv"),sep ='\t',quote=F)
all_samples <- unique(snRNA_obj$sample_ID)
dir.create("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/cell_barcode_tables/")
for (i in 1:length(all_samples)) {
    sample_ID = all_samples[i]
    print(sample_ID)
    sample_barcode_table <- snRNA_obj@meta.data[,c("barcode","cell_type_individual")]
    colnames(sample_barcode_table) <- c("barcode","cell_type_individual")
    write.table(sample_barcode_table, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/",sample_ID,"_cell_types_v7.tsv"),sep='\t',quote=F,row.names=F)
    write.table(sample_barcode_table, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/",sample_ID,"_cell_types_v7.csv"),sep=',',quote=F,row.names=F)
    sample_barcode_table <- snRNA_obj@meta.data[,c("barcode","neoplasm_normal_unknown")]
    colnames(sample_barcode_table) <- c("barcode","neoplasm_normal_unknown")
    write.table(sample_barcode_table, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/",sample_ID,"_neoplasm_normal_unknown_v7.tsv"),sep='\t',quote=F,row.names=F)
    write.table(sample_barcode_table, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/all_snRNA/separate_sct/",sample_ID,"_neoplasm_normal_unknown_v7.csv"),sep=',',quote=F,row.names=F)
}

# Figure 2a and 2c spatial image plots for HT268B1 uses the following code:
# Li wanted me to add a ceiling of one to the plots so that they would be easier to see. That was run in the following manner:
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/
# conda activate seurat5
# Rscript plotting_neoplastic_labels_7_darkblue_ceiling1.R neo_norm_unk_table_v7.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/
# the full section images of the cell segmentation colored by neoplastic vs normal label are from the output of the above script specifically the *_neoplastic_normal_unknown_v4_image_2color.pdf images and the legend comes from *_neoplastic_normal_unknown_v4_image_2color_legend.pdf. by doing it this way I am able to consistently size all of the plots no matter if they are [dimplots or featureplots] or [Imagedimplots or Imagefeatureplots] and then still get the legends. Doing it this way allows the DimPlots and Featureplots to be perfectly square.
# In evaluating that however it was determined that a ceiling of 1 was too low for the UMAP FeaturePlot of EEF1A1 because expression was so high. So I decided to increase the ceiling to 5 for that one as a result.
# the zoom in shown in figure 2c is from where I crop to and create the zoom1 fov below
# the H&E was generated by loading the HEX-SIFT aligned image into Xenium explorer and then exporting the current view with only the H&E image layer active and the cells and transcripts unchecked.
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
options(future.globals.maxSize= +Inf) # need to increase for SP001P1 otherwise it will crash #+Inf is 100% overkill but I don't care right now
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
i = 5 # HT268B1-Th1H3L1U1
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
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
# neoplasm_normal_unknown color codes c('#009E73','#CC79A7','#B4DAF4') ordered by normal, neoplasm, unknown
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic","low_quality/unknown"))
var_key <- "EEF1A1-p-D442H-ALT-G"
ref_key <- "EEF1A1-p-D442H-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
var_tmp <- as.data.frame(counts_df[ ,var_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = var_tmp, col.name = paste0(var_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["EEF1A1_detected"] <- "EEF1A1 p.D442H probe not detected"
obj@meta.data["EEF1A1_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "EEF1A1 p.D442H reference"
obj@meta.data["EEF1A1_detected"][obj@meta.data[paste0(var_key,"_",assay,"_count")] > 0] <- "EEF1A1 p.D442H variant"
obj@meta.data["EEF1A1_detected"][obj@meta.data[paste0(var_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(ref_key,"_",assay,"_count")] ] <- "EEF1A1 p.D442H variant and reference"
obj$EEF1A1_detected <- factor(obj$EEF1A1_detected, levels = c("EEF1A1 p.D442H probe not detected","EEF1A1 p.D442H reference","EEF1A1 p.D442H variant","EEF1A1 p.D442H variant and reference"))
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="EEF1A1_detected", cols = c('darkblue','gold','darkorchid1','firebrick1',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_fov.with.snvs_crop_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="EEF1A1_detected", cols = c('darkblue','gold','darkorchid1',"firebrick1","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_fov.with.snvs_crop_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
zoom1 <- Crop(obj[["fov"]], x = c(4700, 5200), y = c(1900, 2400), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1"]] <- zoom1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="EEF1A1_detected", cols = c('darkblue','gold','darkorchid1',"firebrick1","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_zoom1_crop_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="EEF1A1_detected", cols = c('darkblue','gold','darkorchid1',"firebrick1","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_zoom1_crop_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

coordinates <- Embeddings(obj, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
i = 5 # HT268B1-Th1H3L1U1
sample = input_table[i,"Sample_ID"]
rds_obj = input_table[i,"Xenium_snv_object"]
dir.create(paste0(out_dir,"/",sample))
Manual_cell_type = input_table[i,"Manual_xenium_cell_types_v2"]
neoplastic_labels_v3 <- input_table[i,"Xenium_neoplastic_labels_v3"]
neoplastic_labels_v3_subclone <- input_table[i,"Xenium_neoplastic_labels_v3_subclone"]
unknown_low_quality_labels_v3 <- input_table[i,"Xenium_unknown_low_quality_labels_v3"]
tumor_labels <- unlist(strsplit(neoplastic_labels_v3, ","))
tumor_labels_subclone <- unlist(strsplit(neoplastic_labels_v3_subclone, ","))
unknown_low_quality_labels_v3 <- unlist(strsplit(unknown_low_quality_labels_v3, ","))
feature_csv = input_table[i,"feature_csv"]
feature_list <- unlist(strsplit(feature_csv, ","))
print(sample)
print(feature_list)
obj <- readRDS(rds_obj)
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
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
# neoplasm_normal_unknown color codes c('#009E73','#CC79A7','#B4DAF4') ordered by normal, neoplasm, unknown
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic","low_quality/unknown"))
coordinates <- Embeddings(obj, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
p1 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p2 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p1)
jiterables = length(list_of_plots)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p2
jiterables = length(feature_list)
for (j in 1:jiterables) {
    p3 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], max.cutoff = 2, order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p3
}
print("Feature plot length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_max2.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
# for (j in 1:length(list_of_plots)) {
#     print(list_of_plots[[j]])
# }
# dev.off()
p4 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p5 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p4)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p5
jiterables = length(feature_list)
for (j in 1:jiterables) {
    p6 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], max.cutoff = 2, order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p6
}
print("Feature plot legend length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_max2_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p1 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p2 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p1)
jiterables = length(list_of_plots)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p2
jiterables = length(feature_list)
for (j in 1:jiterables) {
    p3 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], max.cutoff = 5, order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p3
}
print("Feature plot length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_max5.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
# for (j in 1:length(list_of_plots)) {
#     print(list_of_plots[[j]])
# }
# dev.off()
p4 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p5 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p4)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p5
jiterables = length(feature_list)
for (j in 1:jiterables) {
    p6 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], max.cutoff = 5, order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p6
}
print("Feature plot legend length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_max5_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# Figure 2d and 2f spatial image plots for SP001P1-Fp1U1 uses the following code:
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/
# conda activate seurat5
# Rscript plotting_neoplastic_labels_7_darkblue_ceiling1.R neo_norm_unk_table_v7.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/
# the full section images of the cell segmentation colored by neoplastic vs normal label are from the output of the above script specifically the *_neoplastic_normal_unknown_v4_image_2color.pdf images and the legend comes from *_neoplastic_normal_unknown_v4_image_2color_legend.pdf. by doing it this way I am able to consistently size all of the plots no matter if they are [dimplots or featureplots] or [Imagedimplots or Imagefeatureplots] and then still get the legends. Doing it this way allows the DimPlots and Featureplots to be perfectly square.
# the zoom in shown in figure 2c is from where I crop to and create the zoom1.1 fov below.
# there are no normal duct or panin cells in this case
# the H&E was generated by loading the HEX-SIFT aligned image into Xenium explorer and then exporting the current view with only the H&E image layer active and the cells and transcripts unchecked.
# i = 1 # HT227P1-S1H1L1U1 # G12D
# i = 2 # HT242P1-S1H4L4U1 # G12D
i = 19 # SP001P1-Fp1U1 # G12D
# i = 26 # HT270P1-S1H1A1US2_1 # G12D
# i = 27 # SP002C1-Fp1U2 # G12V
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
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
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic","low_quality/unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("normal","normal_duct","PanIN","neoplastic","low_quality/unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
var_key <- "GNAS-p-R201C-ALT-T"
ref_key <- "GNAS-p-R201C-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
var_tmp <- as.data.frame(counts_df[ ,var_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = var_tmp, col.name = paste0(var_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["GNAS_detected"] <- "GNAS p.R201C probe not detected"
obj@meta.data["GNAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "GNAS p.R201C reference"
obj@meta.data["GNAS_detected"][obj@meta.data[paste0(var_key,"_",assay,"_count")] > 0] <- "GNAS p.R201C variant"
obj@meta.data["GNAS_detected"][obj@meta.data[paste0(var_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(ref_key,"_",assay,"_count")] ] <- "GNAS p.R201C variant and reference"
obj$GNAS_detected <- factor(obj$GNAS_detected, levels = c("GNAS p.R201C probe not detected","GNAS p.R201C reference","GNAS p.R201C variant","GNAS p.R201C variant and reference"))
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
zoom1.1 <- Crop(obj[["fov"]], x = c(5500, 6000), y = c(5300, 5800), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="GNAS_detected", cols = c('darkblue','gold','darkorchid1','firebrick1',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_zoom1.1_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="GNAS_detected", cols = c('darkblue','gold','darkorchid1','firebrick1',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_zoom1.1_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="GNAS_detected", cols = c('darkblue','gold','darkorchid1','firebrick1',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_fov.with.snvs_crop_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="GNAS_detected", cols = c('darkblue','gold','darkorchid1',"firebrick1","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_4color_fov.with.snvs_crop_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

#merging the bulk CNV calls from the HTAN and PECGS cases for joint analysis
# needed for figure 3 and related text
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/merge
# conda activate seurat5
library(tidyverse)
## merged.gene_level.from_seg.hg38.calls.tsv
htan <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v3_2025-01-16/merged.gene_level.from_seg.hg38.calls.tsv",sep = "\t", header = 1, row.names = 1)
pecgs <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/PECGS/freeze_v2_2024-07-23/merged.gene_level.from_seg.hg38.calls.tsv",sep = "\t", header = 1, row.names = 1) # there were no new samples from PECGS so I did not re-run this and the results stayed the same from freeze v2
hrow <- rownames(htan)
prow <- rownames(pecgs)
print(hrow[(!(hrow %in% prow))])
# character(0)
print(prow[(!(prow %in% hrow))])
# [1] "OR4F16" "OR4F29"
htan["OR4F16",] <- NA
htan["OR4F29",] <- NA
all_samples <- cbind(htan,pecgs)
dir.create("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/merge/freeze_v3_2025-01-16/")
write.table(all_samples, "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/merge/freeze_v3_2025-01-16/all_samples.gene_level.from_seg.hg38.calls.tsv",sep = "\t",quote=F)
## merged.gene_level.from_seg.hg38.log2ratio.tsv
htan <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v3_2025-01-16/merged.gene_level.from_seg.hg38.log2ratio.tsv",sep = "\t", header = 1, row.names = 1)
pecgs <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/PECGS/freeze_v2_2024-07-23/merged.gene_level.from_seg.hg38.log2ratio.tsv",sep = "\t", header = 1, row.names = 1) # there were no new samples from PECGS so I did not re-run this and the results stayed the same from freeze v2
hrow <- rownames(htan)
prow <- rownames(pecgs)
print(hrow[(!(hrow %in% prow))])
# character(0)
print(prow[(!(prow %in% hrow))])
# [1] "OR4F16" "OR4F29"
htan["OR4F16",] <- NA
htan["OR4F29",] <- NA
all_samples <- cbind(htan,pecgs)
write.table(all_samples, "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/merge/freeze_v3_2025-01-16/all_samples.gene_level.from_seg.hg38.log2ratio.tsv",sep = "\t",quote=F)
## merged.arm_level.from_seg.hg38.calls.tsv
htan <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v3_2025-01-16/merged.arm_level.from_seg.hg38.calls.tsv",sep = "\t", header = 1, row.names = 1)
pecgs <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/PECGS/freeze_v2_2024-07-23/merged.arm_level.from_seg.hg38.calls.tsv",sep = "\t", header = 1, row.names = 1) # there were no new samples from PECGS so I did not re-run this and the results stayed the same from freeze v2
hrow <- rownames(htan)
prow <- rownames(pecgs)
print(hrow[(!(hrow %in% prow))])
# character(0)
print(prow[(!(prow %in% hrow))])
# character(0)
all_samples <- cbind(htan,pecgs)
write.table(all_samples, "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/merge/freeze_v3_2025-01-16/all_samples.arm_level.from_seg.hg38.calls.tsv",sep = "\t",quote=F)
## merged.arm_level.from_seg.hg38.log2ratio.tsv
htan <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/HTAN/freeze_v3_2025-01-16/merged.arm_level.from_seg.hg38.log2ratio.tsv",sep = "\t", header = 1, row.names = 1)
pecgs <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/PECGS/freeze_v2_2024-07-23/merged.arm_level.from_seg.hg38.log2ratio.tsv",sep = "\t", header = 1, row.names = 1) # there were no new samples from PECGS so I did not re-run this and the results stayed the same from freeze v2
hrow <- rownames(htan)
prow <- rownames(pecgs)
print(hrow[(!(hrow %in% prow))])
# character(0)
print(prow[(!(prow %in% hrow))])
# character(0)
all_samples <- cbind(htan,pecgs)
write.table(all_samples, "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/bulk_cnv/merge/freeze_v3_2025-01-16/all_samples.arm_level.from_seg.hg38.log2ratio.tsv",sep = "\t",quote=F)


# For figure 3b and 3f we needed to subset existing merged single cell objects in order to work on just the individual sample HT065B1-
# generating the single cell objects for the HTAN BRCA cases:
# the breast cancer data is split out across 3 objects. 
# 1 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_multiome_snRNA_v10312023_merged_obj_varFeatOnly_annotated.rds
# 2 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_snRNA_merged_obj_celltypeannotation.rds
# 
# 1 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_multiome_snRNA_v10312023_merged_obj_varFeatOnly_annotated_metadata.tsv
# 2 /diskmnt/Projects/HTAN_BRCA_analysis/annotation_test/fernanda/HTAN_BRCA_snRNA_merged_obj_celltypeannotation_metadata.tsv
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
# This was the script for pre-processing and generating the input files to run inferCNV whose output was then plotted on 3b (right).
# conda activate seurat5
#!/usr/bin/env Rscript
library(Seurat)
library(optparse)
require(Matrix)
library(tidyverse)
set.seed(1234)
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="path to seurat rds object",
                metavar="character"),
    make_option(c("-o", "--output"),
                type="character",
                default="./",
                help="parent output folder path (this should be the inputs folder that will be used for inferCNV)",
                metavar="character"),
    make_option(c("-s","--sample_id"),
                type="character",
                default="single_cell_study",
                help="Name of your sample",
                metavar="character"),
    make_option(c("-t","--tumor"),
                type="character",
                default="PDAC",
                help="Name of cancer type - what type of tumor did it come from (e.g. PDAC, UCEC, CRC, GBM, etc.)",
                metavar="character")#,
    # make_option(c("-m","--metadata"),
    #             type="character",
    #             default="/diskmnt/Projects/snATAC_primary/04_celltyped_rds/cell_type_snRNA_merged/v2/CESC/CESC_cell_type.harmonized.cancer.meta.data",
    #             help="spreadsheet containing metadata from multiple seruat objects",
    #             metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#input variables
input=opt$input
sample=opt$sample_id
outPath=opt$output
tumor=opt$tumor
#meta_data_table=opt$metadata
dir.create(paste(outPath,"/raw_counts_matrix/",sep=""), recursive = TRUE)
dir.create(paste(outPath,"/annotations_file/",sample,"/",sep=""), recursive = TRUE)

obj=readRDS(input)
obj <- subset(obj, subset = neoplasm_normal_unknown %in% c("neoplasm","normal"))
DefaultAssay(obj) <- 'RNA'
mat=GetAssayData(object = pbmc, assay = "RNA", slot = "counts")
mat=as.data.frame(mat)
write.table(mat,paste(outPath, "/raw_counts_matrix/", sample,'.raw_counts.tsv',sep=''),sep='\t',quote=FALSE)

tmp_df <- obj@meta.data[,c("original_RNA_barcode","neoplasm_normal_unknown")]
write(paste(sample, "\t", "normal", "\t", tumor, sep=""), file=paste(outPath,"/annotations_file/",'reference_cells.txt',sep=""), append=TRUE)
write.table(tmp_df, paste(outPath,"/annotations_file/",sample,"/",sample,'.Barcode_Annotation.txt',sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

# For Figure 3b and 3c the sample level plots were generated in the following manner:
# remaking the loss of heterozygosity plots using the updated seurat objects and the inferCNV for them
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/freeze_v2_2024-07-22/HT065B1-S1H1A3N1/
# conda activate seurat5
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/freeze_v2_2024-07-22/HT065B1-S1H1A3N1"
sample = "HT065B1-S1H1A3N1"
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1_processed_update_no_doublet_seurat5.0.1.rds") 
infercnv_gene_matrix_results <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/matrices/HT065B1-S1H1A3N1_gene_level.tsv",sep='\t',header=T) # output from inferCNV processing
rownames(infercnv_gene_matrix_results) <- infercnv_gene_matrix_results$index
infercnv_gene_matrix_results$index <- NULL
infercnv_arm_matrix_results <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/matrices/HT065B1-S1H1A3N1_arm_level.tsv", sep='\t', header=T) # output from inferCNV processing
rownames(infercnv_arm_matrix_results) <- infercnv_arm_matrix_results$index
infercnv_arm_matrix_results$index <- NULL
tmp_table <- data.frame(TP53 = infercnv_gene_matrix_results$TP53, row.names = rownames(infercnv_gene_matrix_results))
obj <- AddMetaData(object = obj, metadata = tmp_table, col.name = paste("TP53_copy_ratio"))
tmp_table <- data.frame(chr17p = infercnv_arm_matrix_results$`chr17p`, row.names = rownames(infercnv_arm_matrix_results))
obj <- AddMetaData(object = obj, metadata = tmp_table, col.name = paste("chr17p_copy_ratio"))
tmp_table <- data.frame(chr17q = infercnv_arm_matrix_results$`chr17q`, row.names = rownames(infercnv_arm_matrix_results))
obj <- AddMetaData(object = obj, metadata = tmp_table, col.name = paste("chr17q_copy_ratio"))
tumor_normal_annotation <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/loss_of_heterozygosity/inferCNV/input/annotations_file/HT065B1-S1H1A3N1.Barcode_Annotation.txt", sep='\t', header=F) # input for inferCNV processing
rownames(tumor_normal_annotation) <- tumor_normal_annotation$V1
tumor_normal_annotation$V1 <- NULL
tumor_normal_annotation$V2[tumor_normal_annotation$V2 == "Tumor"] <- "neoplastic"
tumor_normal_annotation$V2[tumor_normal_annotation$V2 == "Normal"] <- "normal"
# color values selected form the okabe and Ito colorblind friendly color palette c("normal","clone_1","clone_2","low_quality/unknown") respectively corresponds to c('#479c76','#c17da5',"#eee461",'#5c5c5c')
tmp_table <- data.frame(tumor_normal = tumor_normal_annotation$V2, row.names = rownames(tumor_normal_annotation))
rownames(tmp_table) <- paste0(sample,"_",rownames(tmp_table))
obj <- AddMetaData(object = obj, metadata = tmp_table, col.name = paste("neoplasm_normal_unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic"))
obj$cell_type_xenium_variant_v1[obj$cell_type_xenium_variant_v1 == "Mono_Macro?"] <- "Mono_Macro"
obj$cell_type_xenium_variant_v1[obj$cell_type_xenium_variant_v1 == "Mast?"] <- "Mast"
table(obj$neoplasm_normal_unknown, obj$TP53_copy_ratio)
#        0.5x   1x
# Normal    0 2240
# Tumor  5021  500
table(obj$neoplasm_normal_unknown, obj$chr17p_copy_ratio)
#        0.5x   1x
# Normal    0 2240
# Tumor  4126 1395
table(obj$neoplasm_normal_unknown, obj$chr17q_copy_ratio)
#          1x
# Normal 2240
# Tumor  5521
coordinates <- Embeddings(obj, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
pdf("DimPlot_HT065B1-S1H1A3N1_TP53_and_chr17pq_copy_ratio.pdf", useDingbats=FALSE, height=5, width = 5)
p1 <- DimPlot(obj, group.by = "TP53_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio) #used the legend from here in 2b (right)
p2 <- DimPlot(obj, group.by = "TP53_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) # used the umap from here in 2b (right)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(obj, group.by = "chr17p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(obj, group.by = "chr17p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(obj, group.by = "chr17q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(obj, group.by = "chr17q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300)) # wrong color values
p1 <- DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(obj, group.by = "cell_type_xenium_variant_v1", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio) # used the umap from here in 2b (left)
p2 <- DimPlot(obj, group.by = "cell_type_xenium_variant_v1", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) # used the umap from here in 2b (left)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
dev.off()
DefaultAssay(obj) <- "SCT"
pdf("FeaturePlot_TP53_expression_HT065B1-S1H1A3N1_snRNA.pdf", useDingbats=FALSE, height=5, width = 5)
p1 <- FeaturePlot(obj, features = "TP53", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- FeaturePlot(obj, features = "TP53", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
dev.off()

# Figure 3d was taken from the output files HT65B1-H1A1A4U1_neoplastic_normal_unknown_v4_image_2color.pdf and the legend from *_neoplastic_normal_unknown_v4_image_2color_legend.pdf generated by running the following commands:
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/
# conda activate seurat5
# mkdir darkblue_ceiling1_v7
# Rscript plotting_neoplastic_labels_7_darkblue_ceiling1.R neo_norm_unk_table_v7.tsv /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/

# For Figure 3e was taken from the output file *_barplots_count_based_counts_p_value_alt_v6.pdf generated by running the following commands:
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/
# conda activate seurat5
# mkdir darkblue_ceiling1_v7
# Rscript plotting_barplots_alt_count_v6.R &> plotting_barplots_alt_count_v6.log

# Figure 3f and 3h analysis/plots were generated below:
# Finding DEGs between the snRNA cancer cells with LOH and those that are still heterozygous for TP53
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(ggrepel)
set.seed(1234)
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1_processed_update_no_doublet_seurat5.0.1.rds")
cell_type_df = read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1_cell_type_xenium_variant_v1.tsv", sep ='\t', header = T)
outdir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
sample_ID = "HT065B1-S1H1A3N1"
dir.create(paste0(outdir,"/",sample_ID))
#infercnv_gene_matrix_results <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/matrices/HT065B1-S1H1A3N1_gene_level.tsv", sep = '\t', header = T)
# already added the inferCNV results and cell types to the metadata in the following table which comes from:
# bash /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/plotting_cnvs_per_cell_type.sh -t plotting_input_table.tsv -o . -g 'ATM,ARID1A,BRIP1,BRCA1,BRCA2,SMAD4,KRAS,TP53,CCND1,LCE1D,CCNE1,MCL1,ITGAE,EOMES,CD200R1,MBP,MYC,EGFR,PDGFRA,ERBB2,TCL1A,CD27'
infercnv_metadata <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/HT065B1-S1H1A3N1_meta.data_with_cnvs.tsv",sep='\t',header=T)
obj@meta.data <- infercnv_metadata
obj$cell_type_tp53_cnv <- paste0(obj$cell_type_xenium_variant_v1, "__", obj$TP53_inferCNV)
# There are two ways to do this.
# One is to run the DEG calculation on the clusters themselves since the CNVs are predominantly grouped by cluster. 
# The other is to run the deg calculation on the tumor cells based on if they are literally called as amplified or deleted.
# We will start with the tumor_cnv method and see what the result is since I view the clusters as a proxy for that in this analysis.
Idents(obj) <- "cell_type_tp53_cnv"
tumor_inferCNV_degs <- FindMarkers(obj,assay = "SCT", slot="data", ident.1= "Tumor__1", ident.2="Tumor__0.5", logfc.threshold = 0, min.pct=0.05, min.diff.pct=0.1, test.use='wilcox')
write.table(tumor_inferCNV_degs,paste0(outdir,sample_ID,"/", sample_ID,"_cell_type_tp53_cnv_Tumor__1_vs_Tumor__0.5_202401001.tsv"),sep="\t", quote=FALSE)
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT65B1-H1A1A4U1/cancer_clone_markers
# Plotting those degs that overlap with Xenium panel geneset 
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT65B1-H1A1A4U1/cancer_clone_markers/"
i = 10 # HT65B1-H1A1A4U1
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
DefaultAssay(obj) <- "SCT" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='SCT') <- 'fov'
DefaultBoundary(obj[["fov"]]) <- "segmentation"
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
# neoplasm_normal_unknown color codes c('#009E73','#CC79A7','#B4DAF4') ordered by normal, neoplasm, unknown
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic","low_quality/unknown"))
# Features from the above DEG analysis comparing heterozygous HT65B1 cancer cells to chr17p del cancer cels that exist in HT065B1 snRNA that overlap the features of the Xenium assay:
degs <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1/HT065B1-S1H1A3N1_cell_type_tp53_cnv_Tumor__1_vs_Tumor__0.5_202401001.tsv", sep='\t',header=T)
degs <- degs[order(degs$avg_log2FC), ]
degs_up_tumor_loh <- degs[(degs$avg_log2FC < 0 ), ]
features_toumor_up_loh <- degs_up_tumor_loh$gene[degs_up_tumor_loh$gene %in% Features(obj)]
pdf(paste0(out_dir,"/",sample,"_deg_Tumor_1_vs_Tumor__0.5_in_Xen_dotplot_up_LOH.pdf"),useDingbats = F, height = 5, width = 15)
p1 <- DotPlot(obj, assay = "SCT", features = features_toumor_up_loh, group.by = "neoplasm_normal_unknown") + RotatedAxis()
print(p1)
dev.off()
pdf(paste0(out_dir,"/",sample,"_deg_Tumor_1_vs_Tumor__0.5_in_Xen_featureplot_up_LOH.pdf"),useDingbats = F,height = 5, width = 5)
for (gene in features_toumor_up_loh) {
    p1 <- FeaturePlot(obj, features = gene, reduction = "umap.30PC", order = F, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
    print(p1)
    p2 <- FeaturePlot(obj, features = gene, reduction = "umap.30PC", order = T, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
    print(p2)
}
dev.off()
degs_down_tumor_down_loh <- degs[(degs$avg_log2FC > 0 ), ]
features_toumor_down_loh <- degs_down_tumor_down_loh$gene[degs_down_tumor_down_loh$gene %in% Features(obj)]
pdf(paste0(out_dir,"/",sample,"_deg_Tumor_1_vs_Tumor__0.5_in_Xen_dotplot_down_LOH.pdf"),useDingbats = F, height = 5, width = 15)
p1 <- DotPlot(obj, assay = "SCT", features = features_toumor_down_loh, group.by = "neoplasm_normal_unknown") + RotatedAxis()
print(p1)
dev.off()
pdf(paste0(out_dir,"/",sample,"_deg_Tumor_1_vs_Tumor__0.5_in_Xen_featureplot_down_LOH.pdf"),useDingbats = F,height = 5, width = 5)
for (gene in features_toumor_down_loh) {
    p1 <- FeaturePlot(obj, features = gene, reduction = "umap.30PC", order = F, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
    print(p1)
    p2 <- FeaturePlot(obj, features = gene, reduction = "umap.30PC", order = T, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
    print(p2)
}
dev.off()
# checking that those degs that are present in the snRNA object are also differentially expressed in the Xenium object
library(future)
options(future.globals.maxSize= +Inf)
plan(multicore, workers = 30)
Idents(obj) <- "cell_type"
xen_degs_cell_type <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(xen_degs_cell_type,paste0(out_dir,"/",sample,"_cell_type_FindAllMarkers_20241008.tsv"),sep="\t", quote=FALSE)
xen_degs_cell_type_normal_duct <- xen_degs_cell_type[(grepl("Normal", xen_degs_cell_type$cluster)),]
xen_degs_cell_type_normal_duct <- xen_degs_cell_type_normal_duct[order(xen_degs_cell_type_normal_duct$pct.1, decreasing=T),]
xen_degs_cell_type_normal_duct <- xen_degs_cell_type_normal_duct[order(xen_degs_cell_type_normal_duct$avg_log2FC, decreasing=T),]
top_20_normal_duct_genes_by_xen <- xen_degs_cell_type_normal_duct[1:20,"gene"]
# plotting the genes that overlap with the snRNA in the snRNA as well:
snRNA <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT065B1-S1H1A3N1_processed_update_no_doublet_seurat5.0.1.rds")
snRNAoutdir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
sample_ID = "HT065B1-S1H1A3N1"
DefaultAssay(snRNA) <- "SCT"
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_FeaturePlots_Tumor_1_vs_Tumor__0.5_up_LOH.pdf"),useDingbats = F,height = 5, width = 5)
for (gene in features_toumor_up_loh) {
    genes_expressed <- Features(obj)
    if (gene %in% genes_expressed) {
        p1 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = F, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p1)
        p2 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = T, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p2)
    }
}
dev.off()
infercnv_metadata <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/plots/HT065B1-S1H1A3N1_meta.data_with_cnvs.tsv",sep='\t',header=T)
snRNA@meta.data <- infercnv_metadata
snRNA$cell_type_tp53_cnv <- paste0(snRNA$cell_type_xenium_variant_v1, "__", snRNA$TP53_inferCNV)
snRNA$normal_tumor_LOH <- "normal"
snRNA$normal_tumor_LOH[snRNA$cell_type_tp53_cnv == "Tumor__0.5"] = "Cancer_TP53_LOH"
snRNA$normal_tumor_LOH[snRNA$cell_type_tp53_cnv == "Tumor__1"] = "Cancer_TP53_heterozygous"
snRNA$normal_tumor_LOH[snRNA$cell_type_xenium_variant_v1 %in% c("Basal_myoepithelial", "Luminal_progenitor")] = "normal_epithelium"
# There are two ways to do this. 
# One is to run the DEG calculation on the clusters themselves since the CNVs are predominantly grouped by cluster. 
# The other is to run the deg calculation on the tumor cells based on if they are literally called as amplified or deleted.
# We will start with the tumor_cnv method and see what the result is since I view the clusters as a proxy for that in this analysis.
Idents(snRNA) <- "normal_tumor_LOH"
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_deg_Tumor_1_vs_Tumor__0.5_in_snRNA_dotplot_up_LOH.pdf"),useDingbats = F, height = 5, width = 15)
p1 <- DotPlot(snRNA, assay = "SCT", features = features_toumor_up_loh, group.by = "normal_tumor_LOH") + RotatedAxis()
print(p1)
dev.off()
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_FeaturePlots_Tumor_1_vs_Tumor__0.5_down_LOH.pdf"),useDingbats = F,height = 5, width = 5)
for (gene in features_toumor_down_loh) {
    genes_expressed <- Features(obj)
    if (gene %in% genes_expressed) {
        p1 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = F, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p1)
        p2 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = T, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p2)
    }
}
dev.off()
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_deg_Tumor_1_vs_Tumor__0.5_in_snRNA_dotplot_down_LOH.pdf"),useDingbats = F, height = 5, width = 15)
p1 <- DotPlot(snRNA, assay = "SCT", features = features_toumor_down_loh, group.by = "normal_tumor_LOH") + RotatedAxis()
print(p1)
dev.off()
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_FeaturePlots_Tumor_1_vs_Tumor__0.5_down_LOH.pdf"),useDingbats = F,height = 5, width = 5)
for (gene in features_toumor_down_loh) {
    genes_expressed <- Features(obj)
    if (gene %in% genes_expressed) {
        p1 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = F, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p1)
        p2 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = T, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p2)
    }
}
dev.off()
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_snRNA_dotplot_of_xen_normal_duct_degs.pdf"),useDingbats = F, height = 5, width = 15)
p1 <- DotPlot(snRNA, assay = "SCT", features = top_20_normal_duct_genes_by_xen, group.by = "normal_tumor_LOH") + RotatedAxis()
print(p1)
dev.off()
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_snRNA_Featureplot_of_xen_normal_duct_degs.pdf"),useDingbats = F,height = 5, width = 5)
for (gene in top_20_normal_duct_genes_by_xen) {
    genes_expressed <- Features(obj)
    if (gene %in% genes_expressed) {
        p1 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = F, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p1)
        p2 <- FeaturePlot(snRNA, features = gene, reduction = "umap.30PC", order = T, label=F, pt.size=0.5, raster=FALSE)# + theme(legend.position = "none"))
        print(p2)
    }
}
dev.off()
# Genes up in LOH that are associated with double-stranded DNA damage repair:
# # BRIP1 - aka BACH1, directly interacts with BRCA1 to facilitate double-stranded break repair PMID: 14983014
# # BRCA1 - orchestrates double-stranded break repair 
# Gene involved in apoptosis regulation:
# # RTKN2 - "Involved in negative regulation of intrinsic apoptotic signaling pathway" - gene cards - PMID: 30389712
# There are several other genes as well that are up in both that seem to be scattered scattered across other signaling/growth pathways
# ELF5
# RTKN2
# ANPEP
# SPDEF
# NPDC1
# SMYD2
# FOXA1
# TBX3
# HES4
# AGR3
# MCF2L
# NTN4
# SMAD4
# BCL2L11
# KRAS
# markers of Luminal breast cancer subtype according to the paper by Michael and Reyka (tumor classified as luminal B based on PAM50 subtyping done by Michael in prior paper):
# FOXA1
# ESR1
# snRNA dot plot 
# Cancer_TP53_LOH
# Cancer_TP53_heterozygous
# normal_epithelium
# normal
snRNA$normal_tumor_LOH <- factor(snRNA$normal_tumor_LOH, levels = c("normal","normal_epithelium","Cancer_TP53_heterozygous","Cancer_TP53_LOH"))
# snRNA dot plot
library(ggh4x)
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_snRNA_dotplot_of_luminal_and_double_strand_break_repair_genes.pdf"),useDingbats = F, height = 1.9, width = 5.4)
p1 <- DotPlot(snRNA, assay = "SCT", features = c("BRIP1","BRCA1","FOXA1","ESR1"), col.min = -0.5, col.max = 1.0, dot.min = 0, dot.scale = 6, group.by = "normal_tumor_LOH") + RotatedAxis()
print(p1)
dev.off()
pdf(paste0(snRNAoutdir,"/",sample_ID,"/",sample_ID,"_snRNA_dotplot_of_luminal_and_double_strand_break_repair_genes_test.pdf"),useDingbats = F, height = 6, width = 6)
p1 <- DotPlot(snRNA, assay = "SCT", features = c("BRIP1","BRCA1","FOXA1","ESR1"), col.min = -0.5, col.max = 1.0, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "normal_tumor_LOH") + RotatedAxis() +
    force_panelsizes(rows = unit(0.88, "in"),
                     cols = unit(0.88, "in"))
print(p1)
dev.off()
obj$normal_tumor_LOH <- "Normal"
obj$normal_tumor_LOH[(obj$cell_type %in% c("Myoepithelial/Normal ducts"))] <- "Normal_epithelium"
obj$normal_tumor_LOH[(obj$cell_type %in% tumor_labels)] <- "Cancer_TP53_LOH"
obj$normal_tumor_LOH <- factor(obj$normal_tumor_LOH, levels = c("Normal","Normal_epithelium","Cancer_TP53_LOH"))
# Xenium dot plot
pdf(paste0(out_dir,"/",sample,"_Xenium_dotplot_of_luminal_and_double_strand_break_repair_genes.pdf"),useDingbats = F, height = 1.8, width = 4.7)
p2 <- DotPlot(obj, assay = "SCT", features = c("BRIP1","BRCA1","FOXA1","ESR1"), col.min = -0.5, col.max = 1.0, dot.min = 0, dot.scale = 6, group.by = "normal_tumor_LOH") + RotatedAxis()
print(p2)
dev.off()
pdf(paste0(out_dir,"/",sample,"_Xenium_dotplot_of_luminal_and_double_strand_break_repair_genes_test.pdf"),useDingbats = F, height = 6, width = 6)
p2 <- DotPlot(obj, assay = "SCT", features = c("BRIP1","BRCA1","FOXA1","ESR1"), col.min = -0.5, col.max = 1.0, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "normal_tumor_LOH") + RotatedAxis() +
    force_panelsizes(rows = unit(0.66, "in"),
                     cols = unit(0.88, "in"))
print(p2)
dev.off()
pdf(paste0(out_dir,"/",sample,"_Xenium_dotplot_of_luminal_and_double_strand_break_repair_genes_test.pdf"),useDingbats = F, height = 6, width = 12)
print(p1|p2)
dev.off()
library(ggpubr)
#use this one in the figure
pdf(paste0(out_dir,"/",sample,"_Xenium_dotplot_of_luminal_and_double_strand_break_repair_genes_testggarrange_3f.pdf"),useDingbats = F, height = 6, width = 12) #this is the dotplot that we used in Figure 3F
print(ggarrange(p1,p2, ncol = 2, nrow = 1))
dev.off()
# the zoomed in spatial cell segmentation in figure 3h was then generated below.
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
zoom2 <- Crop(obj[["fov.with.snvs"]], x = c(1690, 2190), y = c(2000, 2500), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom2"]] <- zoom2 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom2'
DefaultBoundary(obj[["zoom2"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2_crop_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2_crop_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# Figure 3g was taken from the output file *_box_plot_proportion_alt_cells_v6_subclone.pdf generated by running the following commands:
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/
# conda activate seurat5
# mkdir darkblue_ceiling1_v7
# Rscript plotting_boxplots_alt_cells_false_positive_v6_subclone.R

# Figure 4b and Extended data figure 2.
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(Seurat)
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
i = 4 # HT260C1-Th1K1L1U1
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
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
obj$neoplasm_normal_unknown[(obj$cell_type %in% tumor_labels) & (obj$cell_type == "Tumor_1")] <- "clone_1"
obj$neoplasm_normal_unknown[(obj$cell_type %in% tumor_labels) & (obj$cell_type == "Tumor_2")] <- "clone_2"
obj$cell_type[obj$cell_type %in% c('Macrophage_LowCount','Macrophage')] <- 'Macrophage'
obj$cell_type[obj$cell_type %in% c('Fibroblast_LowCount','Fibroblast')] <- 'Fibroblast'
tmp_df <- obj@meta.data[,c("barcode","neoplasm_normal_unknown")]
colnames(tmp_df) <- c("cell_id","neoplasm_normal_unknown_v4_clones")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplasm_normal_unknown_v4_clones.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
tmp_df <- obj@meta.data[,c("barcode","neoplasm_normal_unknown")]
colnames(tmp_df) <- c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplasm_normal_unknown_v4_clones.csv"), row.names=FALSE, sep=",", quote=FALSE)
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","clone_1","clone_2","low_quality/unknown"))
# color values selected form the okabe and Ito colorblind friendly color palette c("normal","clone_1","clone_2","low_quality/unknown") respectively corresponds to c('#009E73','#CC79A7',"#eee461",'#5c5c5c')
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v6_subclone_image_2color_fov.with.snvs_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # this is figure 4b middle
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v6_subclone_image_2color_fov.with.snvs_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # this is the legend for figure 4b middle and right
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
zoom1 <- Crop(obj[["fov"]], x = c(2850, 3350), y = c(1900, 2400), coords = c("plot","tissue")) # 
obj[["zoom1"]] <- zoom1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v6_subclone_image_2color_zoom1_crop_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v6_subclone_image_2color_zoom1_crop_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
zoom2 <- Crop(obj[["fov"]], x = c(3000, 3500), y = c(4750, 5250), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom2"]] <- zoom2 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom2'
DefaultBoundary(obj[["zoom2"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_subclone_image_2color_zoom2_crop_darkblue.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_subclone_image_2color_zoom2_crop_darkblue_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
coordinates <- Embeddings(obj, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
p1 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p2 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p1)
jiterables = length(list_of_plots)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p2
jiterables = length(feature_list)
for (j in 1:jiterables) {
    p3 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p3
}
print("Feature plot length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_subclone.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # this is figure 4b right and 4c
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
# for (j in 1:length(list_of_plots)) {
#     print(list_of_plots[[j]])
# }
# dev.off()
p4 <- rasterize(DimPlot(obj, group.by = "neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p5 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p4)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p5
jiterables = length(feature_list)
for (j in 1:jiterables) {
    p6 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p6
}
print("Feature plot legend length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_subclone_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this is the legend for figure 4b middle and right and the legend for figure 4c
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
library(viridisLite)
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v4_nCount_Xenium_FeaturePlot_legend.pdf"),useDingbats = F,height=5, width = 5) #this is the legend for figure 4b middle and right and the legend for figure 4c
print(rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "nCount_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
print(rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "nCount_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600))
dev.off()

library(epitools)
library(effectsize)
library(scales)
variant_keys <- c("TPD52L1-p-A21---ALT-T","LDHB-p-Q307---ALT-T","APC-p-S1298Ffs--3-ALT-T")
wt_keys <- c("TPD52L1-p-A21---WT","LDHB-p-Q307---WT","APC-p-S1298Ffs--3-WT")
variant_key <- variant_keys[1]
wt_key <- wt_keys[1]
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj$tumor_cell_call <- obj$neoplasm_normal_unknown 
levels(obj$tumor_cell_call) = c("Normal","Cancer_clone_1","Cancer_clone_2","low_quality/unknown")
cell_freq = as.data.frame(table(obj$tumor_cell_call, obj@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
n_cell_clone_1 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_1"),"Freq"])
n_cell_clone_2 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_2"),"Freq"])
n_cell_normal = sum(cell_freq[(cell_freq$Var1 == "Normal"),"Freq"])
data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Normal"))
data_1["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_1["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_1["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_1["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])
#                ALT   WT
# Cancer_clone_2 587 3298
# Normal          47  699
variant_chisq_test_1 <- chisq.test(data_1, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_1)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_1
# X-squared = 41.1, df = NA, p-value = 0.0004998
x <- data_1[,'ALT']
n <- rowSums(data_1)
variant_prop_test_1 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_1)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 41.1, df = 1, p-value = 7.23e-11
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.07067224 1.00000000
# sample estimates:
#     prop 1     prop 2 
# 0.15109395 0.06300268 
oddsratio.out_1 = oddsratio(as.matrix(data_1), method = "fisher")
print("Odds ratio")
print(oddsratio.out_1)
# [1] "Odds ratio"
# Odds ratio |       95% CI
# -------------------------
#     2.65       | [1.95, 3.60]
phi.effectsize_1 <- phi(as.matrix(data_1), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_1)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.09       | [0.07, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Cancer_clone_1"))
data_2["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
data_2["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])
#                ALT    WT
# Cancer_clone_2 587  3298
# Cancer_clone_1  54 11659
variant_chisq_test_2 <- chisq.test(data_2, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_2)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_2
# X-squared = 1588.6, df = NA, p-value = 0.0004998
x <- data_2[,'ALT']
n <- rowSums(data_2)
variant_prop_test_2 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_2)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 1588.6, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.1369766 1.0000000
# sample estimates:
#     prop 1      prop 2 
# 0.151093951 0.004610262 
oddsratio.out_2 = oddsratio(as.matrix(data_2), method = "fisher")
print("Odds ratio")
print(oddsratio.out_2)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     38.43      | [29.00, 50.92]
phi.effectsize_2 <- phi(as.matrix(data_2), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_2)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.32       | [0.31, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_1", "Normal"))
data_3["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_3["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_3["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_3["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
#                ALT    WT
# Cancer_clone_1  54 11659
# Normal          47   699
variant_chisq_test_3 <- chisq.test(data_3, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_3)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_3
# X-squared = 297.39, df = NA, p-value = 0.0004998
x <- data_3[,'ALT']
n <- rowSums(data_3)
variant_prop_test_3 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_3)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 297.39, df = 1, p-value = 1
# alternative hypothesis: greater
# 95 percent confidence interval:
#     -0.07306069  1.00000000
# sample estimates:
#     prop 1      prop 2 
# 0.004610262 0.063002681 
oddsratio.out_3 = oddsratio(as.matrix(data_3), method = "fisher")
print("Odds ratio")
print(oddsratio.out_3)
# [1] "Odds ratio"
# Odds ratio |       95% CI
# -------------------------
#     0.07       | [0.05, 0.10]
phi.effectsize_3 <- phi(as.matrix(data_3), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_3)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.15       | [0.14, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
normal_var_count <- data_1["Normal","ALT"]
normal_ref_count <- data_1["Normal","WT"]
clone_2_var_count <- data_1["Cancer_clone_2","ALT"]
clone_2_ref_count <- data_1["Cancer_clone_2","WT"]
clone_1_var_count <- data_3["Cancer_clone_1","ALT"]
clone_1_ref_count <- data_3["Cancer_clone_1","WT"]

data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      probe_count = c(normal_var_count,normal_ref_count,clone_2_var_count,clone_2_ref_count,clone_1_var_count,clone_1_ref_count),
                      probe_total = rep(c(sum(normal_var_count,normal_ref_count),sum(clone_2_var_count,clone_2_ref_count),sum(clone_1_var_count,clone_1_ref_count)), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
df_max_round = round(max(data_df$probe_total), -2)+100
plot_max = 1.2*df_max_round
cols <- c("Reference" = hue_pal()(length(unique(data_df$Allele_detected)))[2], "Variant" = hue_pal()(length(unique(data_df$Allele_detected)))[1])
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0,1.2*plot_max))+ #, breaks = seq(from = 0, to = 1600, by = 400)) +
    scale_fill_manual(values = cols) +
    # scale_y_continuous(limits = c(1,10000), trans='log10') +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = df_max_round*1.1, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = df_max_round*1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = df_max_round*1.3, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = probe_total),
              color = "black", nudge_y = max(data_df$probe_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Count of ",variant_key," in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_",variant_key,"_barplots_count_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 6)
print(p1)
dev.off()

normal_total = normal_var_count + normal_ref_count
clone_2_total = clone_2_var_count + clone_2_ref_count
clone_1_total = clone_1_var_count + clone_1_ref_count
data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      proportion = c(normal_var_count/normal_total,normal_ref_count/normal_total,clone_2_var_count/clone_2_total,clone_2_ref_count/clone_2_total,clone_1_var_count/clone_1_total,clone_1_ref_count/clone_1_total),
                      prop_total = rep(c(1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    scale_fill_manual(values = cols) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_",variant_key,"_barplots_proportion_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

variant_keys <- c("TPD52L1-p-A21---ALT-T","LDHB-p-Q307---ALT-T","APC-p-S1298Ffs--3-ALT-T")
wt_keys <- c("TPD52L1-p-A21---WT","LDHB-p-Q307---WT","APC-p-S1298Ffs--3-WT")
# LDHB-p-Q307---ALT-T
variant_key <- variant_keys[2]
wt_key <- wt_keys[2]
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj$tumor_cell_call <- obj$neoplasm_normal_unknown 
levels(obj$tumor_cell_call) = c("Normal","Cancer_clone_1","Cancer_clone_2","low_quality/unknown")
cell_freq = as.data.frame(table(obj$tumor_cell_call, obj@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
n_cell_clone_1 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_1"),"Freq"])
n_cell_clone_2 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_2"),"Freq"])
n_cell_normal = sum(cell_freq[(cell_freq$Var1 == "Normal"),"Freq"])
data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Normal"))

data_1["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_1["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_1["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_1["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])
# ALT    WT
# Cancer_clone_2 13424 60274
# Normal           526  8568
variant_chisq_test_1 <- chisq.test(data_1, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_1)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_1
# X-squared = 892.84, df = NA, p-value = 0.0004998
x <- data_1[,'ALT']
n <- rowSums(data_1)
variant_prop_test_1 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_1)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 892.84, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.1196521 1.0000000
# sample estimates:
#     prop 1     prop 2 
# 0.18214877 0.05784033
oddsratio.out_1 = oddsratio(as.matrix(data_1), method = "fisher")
print("Odds ratio")
print(oddsratio.out_1)
# [1] "Odds ratio"
# Odds ratio |       95% CI
# -------------------------
#     3.63       | [3.32, 3.97]
phi.effectsize_1 <- phi(as.matrix(data_1), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_1)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.10       | [0.10, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Cancer_clone_1"))
data_2["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
data_2["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])
# ALT    WT
# Cancer_clone_2 13424 60274
# Cancer_clone_1   830 89372
variant_chisq_test_2 <- chisq.test(data_2, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_2)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_2
# X-squared = 15278, df = NA, p-value = 0.0004998
x <- data_2[,'ALT']
n <- rowSums(data_2)
variant_prop_test_2 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_2)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 15278, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.1705509 1.0000000
# sample estimates:
#     prop 1     prop 2 
# 0.18214877 0.00920157 
oddsratio.out_2 = oddsratio(as.matrix(data_2), method = "fisher")
print("Odds ratio")
print(oddsratio.out_2)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     23.98      | [22.34, 25.74]
phi.effectsize_2 <- phi(as.matrix(data_2), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_2)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.31       | [0.30, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_1", "Normal"))
data_3["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_3["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_3["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_3["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
# ALT    WT
# Cancer_clone_1 830 89372
# Normal         526  8568
variant_chisq_test_3 <- chisq.test(data_3, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_3)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_3
# X-squared = 1450.9, df = NA, p-value = 0.0004998
x <- data_3[,'ALT']
n <- rowSums(data_3)
variant_prop_test_3 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_3)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 1450.9, df = 1, p-value = 1
# alternative hypothesis: greater
# 95 percent confidence interval:
#     -0.05269908  1.00000000
# sample estimates:
#     prop 1     prop 2 
# 0.00920157 0.05784033
oddsratio.out_3 = oddsratio(as.matrix(data_3), method = "fisher")
print("Odds ratio")
print(oddsratio.out_3)
# [1] "Odds ratio"
# Odds ratio |       95% CI
# -------------------------
#     0.15       | [0.14, 0.17]
phi.effectsize_3 <- phi(as.matrix(data_3), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_3)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.12       | [0.12, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
normal_var_count <- data_1["Normal","ALT"]
normal_ref_count <- data_1["Normal","WT"]
clone_2_var_count <- data_1["Cancer_clone_2","ALT"]
clone_2_ref_count <- data_1["Cancer_clone_2","WT"]
clone_1_var_count <- data_3["Cancer_clone_1","ALT"]
clone_1_ref_count <- data_3["Cancer_clone_1","WT"]

data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      probe_count = c(normal_var_count,normal_ref_count,clone_2_var_count,clone_2_ref_count,clone_1_var_count,clone_1_ref_count),
                      probe_total = rep(c(sum(normal_var_count,normal_ref_count),sum(clone_2_var_count,clone_2_ref_count),sum(clone_1_var_count,clone_1_ref_count)), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
df_max_round = round(max(data_df$probe_total), -2)+100
plot_max = 1.2*df_max_round
cols <- c("Reference" = hue_pal()(length(unique(data_df$Allele_detected)))[2], "Variant" = hue_pal()(length(unique(data_df$Allele_detected)))[1])
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0,1.2*plot_max))+ #, breaks = seq(from = 0, to = 1600, by = 400)) +
    scale_fill_manual(values = cols) +
    # scale_y_continuous(limits = c(1,10000), trans='log10') +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = df_max_round*1.1, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = df_max_round*1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = df_max_round*1.3, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = probe_total),
              color = "black", nudge_y = max(data_df$probe_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Count of ",variant_key," in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_",variant_key,"_barplots_count_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 6)
print(p1)
dev.off()

normal_total = normal_var_count + normal_ref_count
clone_2_total = clone_2_var_count + clone_2_ref_count
clone_1_total = clone_1_var_count + clone_1_ref_count
data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      proportion = c(normal_var_count/normal_total,normal_ref_count/normal_total,clone_2_var_count/clone_2_total,clone_2_ref_count/clone_2_total,clone_1_var_count/clone_1_total,clone_1_ref_count/clone_1_total),
                      prop_total = rep(c(1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    scale_fill_manual(values = cols) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_",variant_key,"_barplots_proportion_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()
variant_keys <- c("TPD52L1-p-A21---ALT-T","LDHB-p-Q307---ALT-T","APC-p-S1298Ffs--3-ALT-T")
wt_keys <- c("TPD52L1-p-A21---WT","LDHB-p-Q307---WT","APC-p-S1298Ffs--3-WT")
# APC-p-S1298Ffs--3-ALT-T
variant_key <- variant_keys[3]
wt_key <- wt_keys[3]
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj$tumor_cell_call <- obj$neoplasm_normal_unknown 
levels(obj$tumor_cell_call) = c("Normal","Cancer_clone_1","Cancer_clone_2","low_quality/unknown")
cell_freq = as.data.frame(table(obj$tumor_cell_call, obj@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
n_cell_clone_1 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_1"),"Freq"])
n_cell_clone_2 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_2"),"Freq"])
n_cell_normal = sum(cell_freq[(cell_freq$Var1 == "Normal"),"Freq"])
data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Normal"))
data_1["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_1["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_1["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_1["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])
# ALT  WT
# Cancer_clone_2  84  22
# Normal         101 347
variant_chisq_test_1 <- chisq.test(data_1, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_1)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_1
# X-squared = 123.9, df = NA, p-value = 0.0004998
x <- data_1[,'ALT']
n <- rowSums(data_1)
variant_prop_test_1 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_1)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 123.9, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.4945321 1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.7924528 0.2254464 
oddsratio.out_1 = oddsratio(as.matrix(data_1), method = "fisher")
print("Odds ratio")
print(oddsratio.out_1)
# [1] "Odds ratio"
# Odds ratio |        95% CI
# --------------------------
#     13.12      | [7.81, 22.04]
phi.effectsize_1 <- phi(as.matrix(data_1), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_1)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.47       | [0.40, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Cancer_clone_1"))
data_2["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
data_2["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])
#                ALT WT
# Cancer_clone_2  84 22
# Cancer_clone_1 171 36
variant_chisq_test_2 <- chisq.test(data_2, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_2)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_2
# X-squared = 0.52531, df = NA, p-value = 0.5482
x <- data_2[,'ALT']
n <- rowSums(data_2)
variant_prop_test_2 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_2)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 0.52531, df = 1, p-value = 0.7657
# alternative hypothesis: greater
# 95 percent confidence interval:
#     -0.1115811  1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.7924528 0.8260870
oddsratio.out_2 = oddsratio(as.matrix(data_2), method = "fisher")
print("Odds ratio")
print(oddsratio.out_2)
# [1] "Odds ratio"
# Odds ratio |       95% CI
# -------------------------
#     0.80       | [0.45, 1.45]
phi.effectsize_2 <- phi(as.matrix(data_2), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_2)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.00       | [0.00, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_1", "Normal"))
data_3["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_3["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_3["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_3["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
# ALT  WT
# Cancer_clone_1 171  36
# Normal         101 347
variant_chisq_test_3 <- chisq.test(data_3, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_3)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_3
# X-squared = 210.35, df = NA, p-value = 0.0004998
x <- data_3[,'ALT']
n <- rowSums(data_3)
variant_prop_test_3 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_3)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 210.35, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.5464895 1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.8260870 0.2254464 
oddsratio.out_3 = oddsratio(as.matrix(data_3), method = "fisher")
print("Odds ratio")
print(oddsratio.out_3)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     16.32      | [10.70, 24.89]
phi.effectsize_3 <- phi(as.matrix(data_3), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize_3)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.57       | [0.50, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
normal_var_count <- data_1["Normal","ALT"]
normal_ref_count <- data_1["Normal","WT"]
clone_2_var_count <- data_1["Cancer_clone_2","ALT"]
clone_2_ref_count <- data_1["Cancer_clone_2","WT"]
clone_1_var_count <- data_3["Cancer_clone_1","ALT"]
clone_1_ref_count <- data_3["Cancer_clone_1","WT"]

data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      probe_count = c(normal_var_count,normal_ref_count,clone_2_var_count,clone_2_ref_count,clone_1_var_count,clone_1_ref_count),
                      probe_total = rep(c(sum(normal_var_count,normal_ref_count),sum(clone_2_var_count,clone_2_ref_count),sum(clone_1_var_count,clone_1_ref_count)), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
df_max_round = round(max(data_df$probe_total), -2)+100
plot_max = 1.2*df_max_round
cols <- c("Reference" = hue_pal()(length(unique(data_df$Allele_detected)))[2], "Variant" = hue_pal()(length(unique(data_df$Allele_detected)))[1])
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0,1.2*plot_max))+ #, breaks = seq(from = 0, to = 1600, by = 400)) +
    scale_fill_manual(values = cols) +
    # scale_y_continuous(limits = c(1,10000), trans='log10') +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = df_max_round*1.1, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = df_max_round*1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = df_max_round*1.3, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = probe_total),
              color = "black", nudge_y = max(data_df$probe_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Count of ",variant_key," in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_",variant_key,"_barplots_count_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 6)
print(p1)
dev.off()

normal_total = normal_var_count + normal_ref_count
clone_2_total = clone_2_var_count + clone_2_ref_count
clone_1_total = clone_1_var_count + clone_1_ref_count
data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      proportion = c(normal_var_count/normal_total,normal_ref_count/normal_total,clone_2_var_count/clone_2_total,clone_2_ref_count/clone_2_total,clone_1_var_count/clone_1_total,clone_1_ref_count/clone_1_total),
                      prop_total = rep(c(1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    scale_fill_manual(values = cols) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_",variant_key,"_barplots_proportion_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()
# Figure 4e stacked bar plot after Morph layer calculation from cell_type proportion table output by Andre
# the tab20 colormap from matplotlib (pulled from here: https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2F09509942-4112-474f-9dfb-2698e6a1e4c2%2F26e4e036-74d6-47ba-85f0-6bbd91e51f58%2Ffiles%2Ffunctions%2Fmatplotlib%2Ftab20.m&embed=web)
# it also matches the hex values I get from adobe illustrator when I use the color picker tool
# conda activate seurat5
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(reshape2)
set.seed(1234)
sample_ID = "HT260C1-Th1K1L1U1"
tumor_1_count <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT260C1-Th1K1L1U1/HT260_layer_stacked_bar_plot_tumor_1_source.tsv",sep='\t',header=T,row.names = "cell_type")
tumor_2_count <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT260C1-Th1K1L1U1/HT260_layer_stacked_bar_plot_tumor_2_source.tsv",sep='\t',header=T,row.names = "cell_type")
tab20_colors <- c('#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2','#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5')
tab20_colors_reorder <- c('#9467bd','#1f77b4','#17becf','#ff9896','#d62728','#98df8a','#2ca02c','#e377c2','#ff9896','#ff7f0e','#bcbd22','#2ca02c','#98df8a','#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#17becf','#9edae5','#f7b6d2','#dbdb8d','#c49c94','#7f7f7f','#c7c7c7')
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT260C1-Th1K1L1U1/"
col_names <- colnames(tumor_1_count)
col_names <- gsub("\\.","\\-",col_names)
col_names <- gsub("X","",col_names)
colnames(tumor_1_count) <- col_names
tumor_1_count["total",] = colSums(tumor_1_count)
row_name_vec <- rownames(tumor_1_count)
for (row_name in row_name_vec) {
    tumor_1_count[paste0(row_name,"_proportion"),] <- tumor_1_count[paste0(row_name),]/tumor_1_count["total",]
}
tumor_1_proportion = tumor_1_count[grepl("proportion",rownames(tumor_1_count)),c(1:10)]
tumor_1_proportion = tumor_1_proportion[(!(grepl("total",rownames(tumor_1_proportion)))),]
tumor_1_proportion$cell_type = unlist(strsplit(rownames(tumor_1_proportion), "_proportion"))
rownames(tumor_1_proportion) = tumor_1_proportion$cell_type
layer_labels = col_names[1:10]
tumor_1_proportion_long = melt(tumor_1_proportion, id.vars = c("cell_type"))
colnames(tumor_1_proportion_long) = c("cell_type","layer","proportion")
layer_start_tmp = as.numeric(unlist(strsplit(as.character(tumor_1_proportion_long$layer), "-")))
layer_start_vector <- layer_start_tmp[seq(1, length(layer_start_tmp), 2)]
tumor_1_proportion_long$layer_start <- layer_start_vector
tumor_1_proportion_long <- tumor_1_proportion_long[tumor_1_proportion_long$cell_type %in% c("Macrophage", "Fibroblast", "Endothelial", "Mast", "T cell", "Cholangiocyte","SMC", "Hepatocyte"),]
tumor_1_proportion_long$cell_type <- factor(tumor_1_proportion_long$cell_type, levels = c("Macrophage", "Fibroblast", 
                                                                                          "Endothelial", "Mast", 
                                                                                          "T cell", "Cholangiocyte", 
                                                                                          "SMC", "Hepatocyte"))
write.table(tumor_1_proportion_long, paste0(out_dir,"/","HT260_layer_stacked_bar_plot_tumor_1_proportions_for_plots.tsv"),sep='\t',quote=F)
p1 <- ggplot(tumor_1_proportion_long, aes(x=layer_start, y = proportion, fill = cell_type)) +
    geom_bar(position="fill", stat="identity") + 
    theme_classic() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks=seq(1,300,30),labels=layer_labels) +
    scale_fill_manual( breaks = c("Macrophage", "Fibroblast", 
                                  "Endothelial", "Mast", 
                                  "T cell", "Cholangiocyte", 
                                  "SMC", "Hepatocyte"),
                       values = tab20_colors_reorder[1:8]) +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
pdf(paste0(out_dir,"/",sample_ID,"_Clone_1_source_all_non_cancer_cells_stacked_bar_plot_proportion_layers_color_fix.pdf"),useDingbats = F, height = 5, width = 6)
print(p1)
dev.off()
col_names <- colnames(tumor_2_count)
col_names <- gsub("\\.","\\-",col_names)
col_names <- gsub("X","",col_names)
colnames(tumor_2_count) <- col_names
tumor_2_count["total",] = colSums(tumor_2_count)
row_name_vec <- rownames(tumor_2_count)
for (row_name in row_name_vec) {
    tumor_2_count[paste0(row_name,"_proportion"),] <- tumor_2_count[paste0(row_name),]/tumor_2_count["total",]
}
tumor_2_proportion = tumor_2_count[grepl("proportion",rownames(tumor_2_count)),c(1:10)]
tumor_2_proportion = tumor_2_proportion[(!(grepl("total",rownames(tumor_2_proportion)))),]
tumor_2_proportion$cell_type = unlist(strsplit(rownames(tumor_2_proportion), "_proportion"))
rownames(tumor_2_proportion) = tumor_2_proportion$cell_type
layer_labels = col_names[1:10]
tumor_2_proportion_long = melt(tumor_2_proportion, id.vars = c("cell_type"))
colnames(tumor_2_proportion_long) = c("cell_type","layer","proportion")
layer_start_tmp = as.numeric(unlist(strsplit(as.character(tumor_2_proportion_long$layer), "-")))
layer_start_vector <- layer_start_tmp[seq(1, length(layer_start_tmp), 2)]
tumor_2_proportion_long$layer_start <- layer_start_vector
tumor_2_proportion_long <- tumor_2_proportion_long[tumor_2_proportion_long$cell_type %in% c("Macrophage", "Fibroblast", "Endothelial", "Mast", "T cell", "Cholangiocyte","SMC", "Hepatocyte"),]
tumor_2_proportion_long$cell_type <- factor(tumor_2_proportion_long$cell_type, levels = c("Macrophage", "Fibroblast", 
                                                                                          "Endothelial", "Mast", 
                                                                                          "T cell", "Cholangiocyte", 
                                                                                          "SMC", "Hepatocyte"))
write.table(tumor_2_proportion_long, paste0(out_dir,"/","HT260_layer_stacked_bar_plot_tumor_2_proportions_for_plots.tsv"),sep='\t',quote=F)
p1 <- ggplot(tumor_2_proportion_long, aes(x=layer_start, y = proportion, fill = cell_type)) +
    geom_bar(position="fill", stat="identity") + 
    theme_classic() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks=seq(1,300,30),labels=layer_labels) +
    scale_fill_manual( breaks = c("Macrophage", "Fibroblast", 
                                  "Endothelial", "Mast", 
                                  "T cell", "Cholangiocyte", 
                                  "SMC", "Hepatocyte"),
                       values = tab20_colors_reorder[1:8]) +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
pdf(paste0(out_dir,"/",sample_ID,"_Clone_2_source_all_non_cancer_cells_stacked_bar_plot_proportion_layers_color_fix.pdf"),useDingbats = F, height = 5, width = 6)
print(p1)
dev.off()

# Spatial plots for figure 4f and 4g and the bar plots for 4f
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/
# conda activate seurat5
# Plots of clone1 and clone 2 in C3L-01287-11Us2_1 cancer cells
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
i = 23 # C3L-01287-11Us2_1 # BAP1 p.T93A
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
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
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","neoplastic","low_quality/unknown"))
obj$neoplasm_clone_normal_unknown <- NA
obj$neoplasm_clone_normal_unknown[(obj$cell_type  %in% c("cancer cell"))] <- "Cancer_clone_1"
obj$neoplasm_clone_normal_unknown[(obj$cell_type %in% c("cancer cell region 2"))] <- "Cancer_clone_2"
obj$neoplasm_clone_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_clone_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "Normal epithelial"
obj$neoplasm_clone_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "low_quality/unknown"
obj$neoplasm_clone_normal_unknown <- factor(obj$neoplasm_clone_normal_unknown, levels = c("Normal","Cancer_clone_1","Cancer_clone_2","low_quality/unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_clone_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_labels.csv"),sep=',',quote=F, row.names=F)
# there are normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
# use zoom1 = clone 2
zoom1 <- Crop(obj[["fov"]], x = c(4500, 5000), y = c(5400, 5900), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1"]] <- zoom1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Clone_2" = "#eee461","Clone_1" = "#CC79A7")
DefaultBoundary(obj[["zoom1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15 
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15 
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_image_2color_zoom1_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this is figure 4g (top - rhabdoid clone) middle, right, far right
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_image_2color_zoom1_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

zoom2 <- Crop(obj[["fov"]], x = c(4350, 4850), y = c(3600, 4100), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom2"]] <- zoom2
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom2'
#Simplified_cell_type = c("Normal" = "#009E73","Clone_2" = "#eee461","Clone_1" = "#CC79A7")
DefaultBoundary(obj[["zoom2"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_image_2color_zoom2_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # this is figure 4g (bottom - clear cell clone) middle, right, far right
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_image_2color_zoom2_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()


p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15 
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # This is figure 4f (right)
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "MKI67", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "MYC", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

coordinates <- Embeddings(obj, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
p1 <- rasterize(DimPlot(obj, group.by = "neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p2 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p1)
jiterables = length(list_of_plots)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p2
jiterables = length(feature_list)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p3 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p3
}
print("Feature plot length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
# for (j in 1:length(list_of_plots)) {
#     print(list_of_plots[[j]])
# }
# dev.off()
p4 <- rasterize(DimPlot(obj, group.by = "neoplasm_clone_normal_unknown", cols = c('#009E73','#CC79A7',"#eee461",'#5c5c5c'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p5 <- rasterize(DimPlot(obj, group.by = "cell_type", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
list_of_plots <- list(p4)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p5
jiterables = length(feature_list)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p6 <- rasterize(FeaturePlot(obj, slot="counts", features = feature_list[j], order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+2]] <- p6
}
print("Feature plot legend length")
print(length(list_of_plots))
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_clone_normal_unknown_v3_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots)))
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# creating bar plot comparing allele proportions between samples.
library(epitools)
library(effectsize)
variant_key <- "BAP1-p-T93A-ALT-C"
wt_key <- "BAP1-p-T93A-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj$tumor_cell_call <- obj$neoplasm_clone_normal_unknown 
# obj$tumor_cell_call <- factor(obj$tumor_cell_call, levels = c("Normal","Cancer_clone_1","Cancer_clone_2","low_quality/unknown"))
cell_freq = as.data.frame(table(obj$tumor_cell_call, obj@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
n_cell_clone_1 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_1"),"Freq"])
n_cell_clone_2 = sum(cell_freq[(cell_freq$Var1 == "Cancer_clone_2"),"Freq"])
n_cell_normal = sum(cell_freq[(cell_freq$Var1 == "Normal"),"Freq"])
data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Normal"))

data_1["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_1["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_1["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_1["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])

variant_chisq_test_1 <- chisq.test(data_1, simulate.p.value=TRUE, correct=FALSE)
#variant_chisq_test2 <- chisq.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_1)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_1
# X-squared = 185.56, df = NA, p-value = 0.0004998
x <- data_1[,'ALT']
#print(x)
n <- rowSums(data_1)
#print(n)
variant_prop_test_1 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_1)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 185.56, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.275354 1.000000
# sample estimates:
#     prop 1     prop 2 
# 0.44827586 0.02090137 
oddsratio.out = oddsratio(as.matrix(data_1), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     38.06      | [16.91, 85.67]
phi.effectsize <- phi(as.matrix(data_1), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.34       | [0.30, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_2", "Cancer_clone_1"))

data_2["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])
data_2["Cancer_clone_2","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(wt_key,"_",assay,"_count")])
data_2["Cancer_clone_2","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"),paste0(variant_key,"_",assay,"_count")])

variant_chisq_test_2 <- chisq.test(data_2, simulate.p.value=TRUE, correct=FALSE)
#variant_chisq_test2 <- chisq.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_2)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_2
# X-squared = 195.67, df = NA, p-value = 0.0004998
x <- data_2[,'ALT']
#print(x)
n <- rowSums(data_2)
#print(n)
variant_prop_test_2 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_2)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 195.67, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.277736 1.000000
# sample estimates:
#     prop 1     prop 2 
# 0.44827586 0.01851852
oddsratio.out = oddsratio(as.matrix(data_2), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     43.06      | [18.74, 98.95]
phi.effectsize <- phi(as.matrix(data_2), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.38       | [0.33, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("Cancer_clone_1", "Normal"))

data_3["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_3["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_3["Cancer_clone_1","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(wt_key,"_",assay,"_count")])
data_3["Cancer_clone_1","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_1"),paste0(variant_key,"_",assay,"_count")])

variant_chisq_test_3 <- chisq.test(data_3, simulate.p.value=TRUE, correct=FALSE)
#variant_chisq_test2 <- chisq.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_3)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_3
# X-squared = 0.21004, df = NA, p-value = 0.6882
x <- data_3[,'ALT']
#print(x)
n <- rowSums(data_3)
#print(n)
variant_prop_test_3 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_3)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 0.21004, df = 1, p-value = 0.6766
# alternative hypothesis: greater
# 95 percent confidence interval:
#     -0.01090284  1.00000000
# sample estimates:
#     prop 1     prop 2 
# 0.01851852 0.02090137 
oddsratio.out = oddsratio(as.matrix(data_3), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |       95% CI
# -------------------------
#     0.88       | [0.52, 1.50]
phi.effectsize <- phi(as.matrix(data_3), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.00       | [0.00, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].

normal_var_count <- data_1["Normal","ALT"]
normal_ref_count <- data_1["Normal","WT"]
clone_2_var_count <- data_1["Cancer_clone_2","ALT"]
clone_2_ref_count <- data_1["Cancer_clone_2","WT"]
clone_1_var_count <- data_3["Cancer_clone_1","ALT"]
clone_1_ref_count <- data_3["Cancer_clone_1","WT"]

data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      probe_count = c(normal_var_count,normal_ref_count,clone_2_var_count,clone_2_ref_count,clone_1_var_count,clone_1_ref_count),
                      probe_total = rep(c(sum(normal_var_count,normal_ref_count),sum(clone_2_var_count,clone_2_ref_count),sum(clone_1_var_count,clone_1_ref_count)), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0,2000), breaks = seq(from = 0, to = 1600, by = 400)) +
    # scale_y_continuous(limits = c(1,10000), trans='log10') +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1700, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1800, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1900, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = probe_total),
              color = "black", nudge_y = max(data_df$probe_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Count of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_count_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

normal_total = normal_var_count + normal_ref_count
clone_2_total = clone_2_var_count + clone_2_ref_count
clone_1_total = clone_1_var_count + clone_1_ref_count
data_df <- data.frame(cell_type = rep(c("Normal cell", "Cancer clone 2 cell", "Cancer clone 1 cell"), each = 2),
                      proportion = c(normal_var_count/normal_total,normal_ref_count/normal_total,clone_2_var_count/clone_2_total,clone_2_ref_count/clone_2_total,clone_1_var_count/clone_1_total,clone_1_ref_count/clone_1_total),
                      prop_total = rep(c(1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell","Cancer clone 2 cell","Cancer clone 1 cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 2 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Cancer clone 2 cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "Cancer clone 1 cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_proportion_p_value_alt_subclones.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

# this shows that cluster 8 in the Seurat object is the rhabdoid cancer clone from this paper: PMID:36563681 see figure 3E and the associated text. This is the same case.
# clusters 0,1,2 are the clear cell renal cell carcinoma.
# based on this result I can then look for DEGs in the single nuc seurat object and see if there is any consistency in the xenium object.
snRNA <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/CPT0079400003_processed_update_no_doublet_seurat5.0.1.rds")
snRNA_sample <- "CPT0079400003"
# infercnv_gene_matrix_results <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/matrices/CPT0079400003_gene_level.tsv",sep='\t',header=T)
# rownames(infercnv_gene_matrix_results) <- infercnv_gene_matrix_results$index
# infercnv_gene_matrix_results$index <- NULL
infercnv_arm_matrix_results <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/outputs/matrices/matrix_cutoffs/CPT0079400003_arm_level.tsv", sep='\t', header=T)
rownames(infercnv_arm_matrix_results) <- infercnv_arm_matrix_results$index
infercnv_arm_matrix_results$index <- NULL
# tmp_table <- data.frame(MYC = infercnv_gene_matrix_results$`MYC`, row.names = rownames(infercnv_gene_matrix_results)) # MYC is not present in the snRNA object inferCNV output (not enough coverage)
# snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("MYC_copy_ratio"))
tmp_table <- data.frame(chr3p = infercnv_arm_matrix_results$`chr3p`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr3p_copy_ratio"))
tmp_table <- data.frame(chr3q = infercnv_arm_matrix_results$`chr3q`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr3q_copy_ratio"))
tmp_table <- data.frame(chr8p = infercnv_arm_matrix_results$`chr8p`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr8p_copy_ratio"))
tmp_table <- data.frame(chr8q = infercnv_arm_matrix_results$`chr8q`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr8q_copy_ratio"))
tmp_table <- data.frame(chr2p = infercnv_arm_matrix_results$`chr2p`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr2p_copy_ratio"))
tmp_table <- data.frame(chr2q = infercnv_arm_matrix_results$`chr2q`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr2q_copy_ratio"))
tmp_table <- data.frame(chr5p = infercnv_arm_matrix_results$`chr5p`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr5p_copy_ratio"))
tmp_table <- data.frame(chr5q = infercnv_arm_matrix_results$`chr5q`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr5q_copy_ratio"))
tmp_table <- data.frame(chr14q = infercnv_arm_matrix_results$`chr14q`, row.names = rownames(infercnv_arm_matrix_results))
snRNA <- AddMetaData(object = snRNA, metadata = tmp_table, col.name = paste("chr14q_copy_ratio"))
tumor_normal_annotation <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/inferCNV/inputs/annotations_file/CPT0079400003/CPT0079400003.Barcode_Annotation.txt", sep='\t', header=F)
rownames(tumor_normal_annotation) <- tumor_normal_annotation$V1
tumor_normal_annotation$V1 <- NULL
tumor_normal_annotation$V2[tumor_normal_annotation$V2 == "neoplasm"] <- "neoplastic"
tumor_normal_annotation$V2[tumor_normal_annotation$V2 == "Normal"] <- "normal"
# color values selected form the okabe and Ito colorblind friendly color palette c("normal","clone_1","clone_2","low_quality/unknown") respectively corresponds to c('#479c76','#c17da5',"#eee461",'#5c5c5c')
snRNA <- AddMetaData(object = snRNA, metadata = tumor_normal_annotation, col.name = paste("infercnv_neoplasm_normal_unknown"))
snRNA$infercnv_neoplasm_normal_unknown <- factor(snRNA$infercnv_neoplasm_normal_unknown, levels = c("normal","neoplastic"))
# table(snRNA$infercnv_neoplasm_normal_unknown, snRNA$MYC_copy_ratio)
table(snRNA$infercnv_neoplasm_normal_unknown, snRNA$chr3p_copy_ratio)
#             0.5 0.529239766081871    1
# neoplasm 3787               145    0
# normal      0                 0 2891
table(snRNA$infercnv_neoplasm_normal_unknown, snRNA$chr3q_copy_ratio)
#          1x
# Normal 2240
# Tumor  5521

snRNA$cell_type_xenium_variant_v2 <- snRNA$cell_type_xenium_variant_v1
snRNA$cell_type_xenium_variant_v2[(snRNA$seurat_clusters == 8) & (snRNA$cell_type_xenium_variant_v2 == "ccRCC cancer cell")] <- "rhabdoid cancer cell"
unique(snRNA$cell_type_xenium_variant_v2)
# [1] "ccRCC cancer cell"            "Loop of Henle"               
# [3] "CD4 T"                        "imyCAF"                      
# [5] "Macrophages_10"               "EX CD8 T"                    
# [7] "Pericytes"                    "CD8 T"                       
# [9] "Macrophages_5"                "NK"                          
# [11] "Macrophages_2"                "Endothelial"                 
# [13] "Macrophages_0"                "rhabdoid cancer cell"        
# [15] "Macrophages_8"                "Macrophages_6"               
# [17] "Treg"                         "CD4CD8 T"                    
# [19] "Macrophages_1"                "Macrophages_3"               
# [21] "mregDC"                       "CAF"                         
# [23] "B cells"                      "Macrophages_7"               
# [25] "MAIT"                         "Macrophages_12"              
# [27] "Plasma cells"                 "Macrophages_9"               
# [29] "Macrophages_4"                "T follicular helpers"        
# [31] "Conventional dendritic cells" "Interferon response T"       
# [33] "Macrophages_11"               "Granulocytes"                
# [35] "Podocyte"                     "Distal tubule"               
# [37] "Proximal tubule"             
snRNA$plotting_reduced_cell_type_xenium_variant_v2 <- snRNA$cell_type_xenium_variant_v2
snRNA$plotting_reduced_cell_type_xenium_variant_v2[ snRNA$cell_type_xenium_variant_v2 %in% c('Macrophages_10','Macrophages_5','Macrophages_2','Macrophages_0','Macrophages_8','Macrophages_6','Macrophages_1','Macrophages_3','Macrophages_7','Macrophages_12','Macrophages_9','Macrophages_4','Macrophages_11')] <- "Macrophages"
snRNA$plotting_reduced_cell_type_xenium_variant_v2[ snRNA$cell_type_xenium_variant_v2 %in% c('CD4 T','CD8 T','EX CD8 T','MAIT','T follicular helpers','CD4CD8 T','Interferon response T','Treg')] <- "T cells"
snRNA$plotting_reduced_cell_type_xenium_variant_v2[ snRNA$cell_type_xenium_variant_v2 %in% c('Conventional dendritic cells','mregDC')] <- "dendritic cells"
snRNA$plotting_reduced_cell_type_xenium_variant_v2[ snRNA$cell_type_xenium_variant_v2 %in% c('Distal tubule','Loop of Henle','Podocyte','Proximal tubule')] <- "normal epithelium"
snRNA$plotting_reduced_cell_type_xenium_variant_v2[ snRNA$cell_type_xenium_variant_v2 %in% c('CAF','imyCAF')] <- "CAF"
Idents(snRNA) <- "plotting_reduced_cell_type_xenium_variant_v2"
snRNA_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/"
tmp_df <- snRNA@meta.data[,c("original_RNA_barcode","plotting_reduced_cell_type_xenium_variant_v2")]
write.table(tmp_df, paste0(snRNA_dir,"/",snRNA_sample,"_plotting_reduced_cell_type_xenium_variant_v2.tsv"),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

coordinates <- Embeddings(snRNA, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
pdf(paste0(snRNA_dir,"/",snRNA_sample,"/20250111/DimPlot_",snRNA_sample,"_chr_arm_copy_ratios_from_PMID36563681_figure3E.pdf"), useDingbats=FALSE, height=5, width = 5)
# p1 <- DimPlot(snRNA, group.by = "MYC_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
# p2 <- DimPlot(snRNA, group.by = "MYC_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
# p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
# print(rasterize(p1, layers='Point', dpi=300))
# print(rasterize(p2, layers='Point', dpi=300))

p1 <- DimPlot(snRNA, group.by = "chr14q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr14q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))

p1 <- DimPlot(snRNA, group.by = "chr3p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr3p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(snRNA, group.by = "chr3q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr3q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))

p1 <- DimPlot(snRNA, group.by = "chr8p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr8p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(snRNA, group.by = "chr8q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr8q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))

p1 <- DimPlot(snRNA, group.by = "chr2p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr2p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(snRNA, group.by = "chr2q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr2q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))

p1 <- DimPlot(snRNA, group.by = "chr5p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr5p_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(snRNA, group.by = "chr5q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "chr5q_copy_ratio", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))

p1 <- DimPlot(snRNA, group.by = "infercnv_neoplasm_normal_unknown", cols = c('#009E73','#CC79A7'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "infercnv_neoplasm_normal_unknown", cols = c('#009E73','#CC79A7'), reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- DimPlot(snRNA, group.by = "plotting_reduced_cell_type_xenium_variant_v2", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- DimPlot(snRNA, group.by = "plotting_reduced_cell_type_xenium_variant_v2", reduction = "umap.30PC", pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
dev.off()
DefaultAssay(snRNA) <- "SCT"
pdf(paste0(snRNA_dir,"/",snRNA_sample,"/20250111/FeaturePlot_expression_",snRNA_sample,"_snRNA_genes_from_PMID36563681_figure3E.pdf"), useDingbats=FALSE, height=5, width = 5)
p1 <- FeaturePlot(snRNA, features = "CA9", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- FeaturePlot(snRNA, features = "CA9", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- FeaturePlot(snRNA, features = "MYC", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- FeaturePlot(snRNA, features = "MYC", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- FeaturePlot(snRNA, features = "MKI67", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- FeaturePlot(snRNA, features = "MKI67", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- FeaturePlot(snRNA, features = "KIF2A", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- FeaturePlot(snRNA, features = "KIF2A", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
p1 <- FeaturePlot(snRNA, features = "MET", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + theme(legend.position="none") + coord_fixed(ratio=xy_ratio)
p2 <- FeaturePlot(snRNA, features = "MET", reduction = "umap.30PC", order = T, label=F, pt.size=0.3, raster=FALSE) + coord_fixed(ratio=xy_ratio)
p1 <- p1 + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
print(rasterize(p1, layers='Point', dpi=300))
print(rasterize(p2, layers='Point', dpi=300))
dev.off()
write.table(snRNA@meta.data, paste0(snRNA_dir,"/",snRNA_sample,"/20250111/",snRNA_sample,"_metadata_20241107.tsv"),sep='\t',quote=F)
library(ggh4x)
library(ggpubr)
snRNA$dotplot_reduced_cell_type <- snRNA$plotting_reduced_cell_type_xenium_variant_v2
snRNA$dotplot_reduced_cell_type[!(snRNA$dotplot_reduced_cell_type %in% c("ccRCC cancer cell","rhabdoid cancer cell"))] <- "Normal"
snRNA$dotplot_reduced_cell_type[snRNA$dotplot_reduced_cell_type %in% c("ccRCC cancer cell")] <- "ccRCC cancer"
snRNA$dotplot_reduced_cell_type[snRNA$dotplot_reduced_cell_type %in% c("rhabdoid cancer cell")] <- "Rhabdoid cancer"
snRNA$dotplot_reduced_cell_type <- factor(snRNA$dotplot_reduced_cell_type, levels = c("Normal", "ccRCC cancer", "Rhabdoid cancer"))
pdf(paste0(snRNA_dir,"/",snRNA_sample,"/20250111/DotPlot_expression_",snRNA_sample,"_snRNA_genes_from_Ilya_pathway_analysis_figure4H.pdf"), useDingbats=FALSE, height=6, width = 6)
p1 <- DotPlot(snRNA, assay = "SCT", features = c("BRCA1","BRIP1","MYC","CENPF","TOP2A","UBE2C","MKI67","BBOX1","GATM","NAT8"), col.min = -0.5, col.max = 2, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "dotplot_reduced_cell_type") + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red") + RotatedAxis() +
    force_panelsizes(rows = unit(0.6, "in"),
                     cols = unit(1.6, "in"))
print(p1)
dev.off()
# obj$normal_tumor_LOH <- "Normal"
# obj$normal_tumor_LOH[(obj$cell_type %in% c("Myoepithelial/Normal ducts"))] <- "Normal_epithelium"
# obj$normal_tumor_LOH[(obj$cell_type %in% tumor_labels)] <- "Cancer_TP53_LOH"
# obj$normal_tumor_LOH <- factor(obj$normal_tumor_LOH, levels = c("Normal","Normal_epithelium","Cancer_TP53_LOH"))
# Xenium dot plot
pdf(paste0(out_dir,"/",sample,"/",sample,"_Xenium_dotplot_genes_from_Ilya_pathway_analysis_figure4H.pdf"),useDingbats = F, height = 6, width = 6)
p2 <- DotPlot(obj, assay = "SCT", features = c("BRCA1","BRIP1","MYC","CENPF","TOP2A","UBE2C","MKI67","BBOX1","GATM","NAT8"), col.min = -0.5, col.max = 2, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "tumor_cell_call") + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red") + RotatedAxis() +
    force_panelsizes(rows = unit(0.6, "in"),
                     cols = unit(1.6, "in"))
print(p2)
dev.off()
#use this one in the figure
pdf(paste0(out_dir,"/",sample,"/",sample,"_Xenium_dotplot_genes_from_Ilya_pathway_analysis_ggarrange_with_snRNA_4h.pdf"),useDingbats = F, height = 6, width = 12)  #this is the dotplots that we used in Figure 4h
print(ggarrange(p1,p2, ncol = 2, nrow = 1))
dev.off()

# Figure 5b
# conda activate seurat5
# plotting the zoomed in images for figure 5 related supplement
# Plotting features for the pre-cancer progression part zoomed in on regions of normal, PanIN, and PDAC
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
# this is the panin from case HT061P1 piece HT061P1-S1P1A1
i = 29 # HT061P1-S1P1A1L1U1
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# This sample was generated using a custom panel where a few of the genes were 
# not present on the panel that were used on  some of the other samples (HT270P1 
# and SP002C1). As a result, during the merged
# analysis, I had to filter the cell by gene matrix of this object to only include
# the genes that were present in all assays. Because some genes were as a result
# dropped from this objects count matrix certain cells then had zero counts. These
# zero count cells had to be filtered out of the object prior to normalization and
# were thus not present during the cell typing of the object. There were 39 cells 
# where this was the case for HT270P1-S1H1A1US2_1. Because this all resulted from these
# 39 cells having low counts they will now be labeled as "LowCount" for future analysis.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type == "Tumor" )] <- "PDAC"
obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_normal_unknown[(obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive"))] <- "Normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Normal_duct","PDAC","PanIN","LowCount/Unknown"))
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS_p.G12_probe_not_detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS_p.G12_reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS_p.G12V_variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS_p.G12D_variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS_p.G12D_p.G12V_variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS_p.G12_probe_not_detected","KRAS_p.G12_reference","KRAS_p.G12V_variant","KRAS_p.G12D_variant","KRAS_p.G12D_p.G12V_variant"))
zoom1.1 <- Crop(obj[["fov"]], x = c(590, 990), y = c(4140,4540), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#eee461',"#56B4E9",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b middle for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#eee461',"#56B4E9",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b middle for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# The fov.with.snvs is in supplement.
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73',"#56B4E9",'#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="KRAS_detected", cols = c('darkblue','yellow','green',"red","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b middle for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73',"#56B4E9",'#CC79A7','#eee461','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="KRAS_detected", cols = c('darkblue','yellow','green',"red","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b middle for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

# this is the tumor and normal duct from case HT061P1 piece HT061P1-S1P1A1
i = 30 # HT061P1-S1P1A1L4U1
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# This sample was generated using a custom panel where a few of the genes were 
# not present on the panel that were used on  some of the other samples (HT270P1 
# and SP002C1). As a result, during the merged
# analysis, I had to filter the cell by gene matrix of this object to only include
# the genes that were present in all assays. Because some genes were as a result
# dropped from this objects count matrix certain cells then had zero counts. These
# zero count cells had to be filtered out of the object prior to normalization and
# were thus not present during the cell typing of the object. There were 39 cells 
# where this was the case for HT270P1-S1H1A1US2_1. Because this all resulted from these
# 39 cells having low counts they will now be labeled as "LowCount" for future analysis.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type == "Tumor" )] <- "PDAC"
obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_normal_unknown[(obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive"))] <- "Normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Normal_duct","PDAC","PanIN","LowCount/Unknown"))
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS_p.G12_probe_not_detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS_p.G12_reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS_p.G12V_variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS_p.G12D_variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS_p.G12D_p.G12V_variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS_p.G12_probe_not_detected","KRAS_p.G12_reference","KRAS_p.G12V_variant","KRAS_p.G12D_variant","KRAS_p.G12D_p.G12V_variant"))
zoom1.1 <- Crop(obj[["fov"]], x = c(1350,1750), y = c(2600, 3000), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#56B4E9','#5c5c5c',"#CC79A7",'#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b left
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#56B4E9','#5c5c5c',"#CC79A7",'#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b left for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
zoom2.1 <- Crop(obj[["fov"]], x = c(3950, 4350), y = c(2750,3150), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom2.1"]] <- zoom2.1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom2.1"]]) <- "segmentation"
# The colors are based on the order in which the factor occurs in the metadata. 
# Additionally some factor levels are not present in the zoom1.2 fov. 
# I use a differently ordered color vector to ensure that colors are consistent across plots. 
# list("normal" = '#009E73',"normal_duct" = "#56B4E9", "PDAC" = "#CC79A7", "PanIN" = "#eee461", "low_quality/unknown" = "#5c5c5c")
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#5c5c5c",'#56B4E9','#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2.1_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b right
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7',"#5c5c5c",'#56B4E9','#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2.1_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b right legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# The fov.with.snvs is in supplement.
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73',"#56B4E9",'#CC79A7','#5c5c5c','#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="KRAS_detected", cols = c('darkblue','yellow','green',"red","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b middle for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_normal_unknown", cols = c('#009E73',"#56B4E9",'#CC79A7','#5c5c5c','#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="KRAS_detected", cols = c('darkblue','yellow','green',"red","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p16
}
DefaultAssay(obj) <- "SCT"
p16 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
p16 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p16
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1_legend.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in figure 5b middle for the legend
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

i = 29 # HT061P1-S1P1A1L1U1
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type == "Tumor" )] <- "PDAC"
obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_normal_unknown[(obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive"))] <- "Normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Normal_duct","PDAC","PanIN","LowCount/Unknown"))

i = 30 # HT061P1-S1P1A1L4U1
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
obj_2 <- readRDS(rds_obj)
DefaultAssay(obj_2) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj_2) <- "seurat_clusters"
DefaultFOV(obj_2, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj_2[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj_2$barcode <- rownames(obj_2@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj_2 <- AddMetaData(object = obj_2, metadata = cell_types, col.name = "cell_type")
obj_2$cell_type[is.na(obj_2$cell_type)] <- "LowCount"
obj_2$neoplasm_normal_unknown <- NA
obj_2$neoplasm_normal_unknown[(obj_2$cell_type == "Tumor" )] <- "PDAC"
obj_2$neoplasm_normal_unknown[(obj_2$cell_type == "PanIN" )] <- "PanIN"
obj_2$neoplasm_normal_unknown[(!(obj_2$cell_type %in% tumor_labels))] <- "Normal"
obj_2$neoplasm_normal_unknown[(obj_2$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive"))] <- "Normal_duct"
obj_2$neoplasm_normal_unknown[(obj_2$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj_2$neoplasm_normal_unknown <- factor(obj_2$neoplasm_normal_unknown, levels = c("Normal","Normal_duct","PDAC","PanIN","LowCount/Unknown"))
# HT061P1-S1P1A1
sample = "HT061P1-S1P1A1"
library(epitools)
library(effectsize)
variant_key <- "KRAS-p-G12D-ALT-T"
wt_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj$tumor_cell_call <- obj$neoplasm_normal_unknown 
# obj$tumor_cell_call <- factor(obj$tumor_cell_call, levels = c("Normal","PanIN","Cancer_clone_2","low_quality/unknown"))
cell_freq = as.data.frame(table(obj$tumor_cell_call, obj@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
assay_class = class(obj_2[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj_2, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj_2 <- AddMetaData(object = obj_2, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj_2 <- AddMetaData(object = obj_2, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj_2@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj_2@meta.data[paste0(variant_key,"_",assay,"_call")][obj_2@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj_2@meta.data[paste0(variant_key,"_",assay,"_call")][obj_2@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj_2$tumor_cell_call <- obj_2$neoplasm_normal_unknown 
# obj$tumor_cell_call <- factor(obj$tumor_cell_call, levels = c("Normal","PanIN","Cancer_clone_2","low_quality/unknown"))
cell_freq_2 = as.data.frame(table(obj_2$tumor_cell_call, obj_2@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
for (var1 in unique(cell_freq$Var1)) {
    for (var2 in unique(cell_freq_2$Var2)) {
        cell_freq$Freq[(cell_freq$Var1 == var1) & (cell_freq$Var2 == var2)] = cell_freq_2$Freq[(cell_freq_2$Var1 == var1) & (cell_freq_2$Var2 == var2)] + cell_freq$Freq[(cell_freq$Var1 == var1) & (cell_freq$Var2 == var2)]
    }
}
n_cell_panin = sum(cell_freq[(cell_freq$Var1 == "PanIN"),"Freq"])
n_cell_pdac = sum(cell_freq[(cell_freq$Var1 == "PDAC"),"Freq"])
n_cell_normal = sum(cell_freq[(cell_freq$Var1 == "Normal"),"Freq"])
n_cell_normal_duct = sum(cell_freq[(cell_freq$Var1 == "Normal_duct"),"Freq"])
data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PDAC", "Normal"))
cell_freq$Freq[cell_freq$Var1 == "Normal" & cell_freq$Var2 == ALT]
data_1["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"), paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal"), paste0(wt_key,"_",assay,"_count")]) 
data_1["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"), paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal"), paste0(variant_key,"_",assay,"_count")])
data_1["PDAC","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PDAC"), paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PDAC"), paste0(wt_key,"_",assay,"_count")])
data_1["PDAC","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PDAC"), paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PDAC"), paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_1 <- chisq.test(data_1, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_1)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_1
# X-squared = 985.62, df = NA, p-value = 0.0004998
x <- data_1[,'ALT']
n <- rowSums(data_1)
variant_prop_test_1 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_1)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 985.62, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.5102803 1.0000000
# sample estimates:
#     prop 1     prop 2 
# 0.62596600 0.09085747 
oddsratio.out = oddsratio(as.matrix(data_1), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |        95% CI
# --------------------------
#     16.75      | [13.74, 20.41]
phi.effectsize <- phi(as.matrix(data_1), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.57       | [0.54, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PDAC", "PanIN"))
data_2["PanIN","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")])
data_2["PanIN","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")])
data_2["PDAC","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PDAC"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PDAC"),paste0(wt_key,"_",assay,"_count")])
data_2["PDAC","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PDAC"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PDAC"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_2 <- chisq.test(data_2, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_2)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_2
# X-squared = 27.62, df = NA, p-value = 0.0004998
x <- data_2[,'ALT']
n <- rowSums(data_2)
variant_prop_test_2 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_2)
#[1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 27.62, df = 1, p-value = 7.381e-08
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.1174732 1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.6259660 0.4542125
oddsratio.out = oddsratio(as.matrix(data_2), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |        95% CI
# --------------------------
#     2.01       | [1.55, 2.62]
phi.effectsize <- phi(as.matrix(data_2), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.13       | [0.09, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PanIN", "Normal"))
data_3["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_3["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_3["PanIN","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")])
data_3["PanIN","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_3 <- chisq.test(data_3, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_3)
# # [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_3
# X-squared = 259.76, df = NA, p-value = 0.0004998
x <- data_3[,'ALT']
n <- rowSums(data_3)
variant_prop_test_3 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_3)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 259.76, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.3125245 1.0000000
# sample estimates:
#     prop 1     prop 2 
# 0.45421245 0.09085747
oddsratio.out = oddsratio(as.matrix(data_3), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |          95% CI
# ----------------------------
#     8.33       | [6.24, 11.11]
phi.effectsize <- phi(as.matrix(data_3), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.36       | [0.32, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_4 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PanIN", "Normal duct"))
data_4["Normal duct","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_duct"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal_duct"),paste0(wt_key,"_",assay,"_count")])
data_4["Normal duct","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_duct"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal_duct"),paste0(variant_key,"_",assay,"_count")])
data_4["PanIN","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")])
data_4["PanIN","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_4 <- chisq.test(data_4, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_4)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_4
# X-squared = 133.74, df = NA, p-value = 0.0004998
x <- data_4[,'ALT']
n <- rowSums(data_4)
variant_prop_test_4 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_4)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 133.74, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.3791502 1.0000000
# sample estimates:
#     prop 1     prop 2 
# 0.45421245 0.02316602 
oddsratio.out = oddsratio(as.matrix(data_4), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     35.09      | [15.09, 81.60]
phi.effectsize <- phi(as.matrix(data_4), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.50       | [0.43, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_5 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PDAC", "Normal duct"))
data_5["Normal duct","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_duct"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal_duct"),paste0(wt_key,"_",assay,"_count")])
data_5["Normal duct","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_duct"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "Normal_duct"),paste0(variant_key,"_",assay,"_count")])
data_5["PDAC","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PDAC"),paste0(wt_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PDAC"),paste0(wt_key,"_",assay,"_count")])
data_5["PDAC","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PDAC"),paste0(variant_key,"_",assay,"_count")]) + sum(obj_2@meta.data[(obj_2@meta.data$tumor_cell_call == "PDAC"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_5 <- chisq.test(data_5, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_5)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_5
# X-squared = 314.48, df = NA, p-value = 0.0004998
x <- data_5[,'ALT']
n <- rowSums(data_5)
variant_prop_test_5 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_5)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 314.48, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.575857 1.000000
# sample estimates:
#     prop 1     prop 2 
# 0.62596600 0.02316602 
oddsratio.out = oddsratio(as.matrix(data_5), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |           95% CI
# -----------------------------
#     70.57      | [31.16, 159.81]
phi.effectsize <- phi(as.matrix(data_5), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.45       | [0.41, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
normal_var_count <- data_1["Normal","ALT"]
normal_ref_count <- data_1["Normal","WT"]
pdac_var_count <- data_1["PDAC","ALT"]
pdac_ref_count <- data_1["PDAC","WT"]
panin_var_count <- data_3["PanIN","ALT"]
panin_ref_count <- data_3["PanIN","WT"]
normal_duct_var_count <- data_4["Normal duct","ALT"]
normal_duct_ref_count <- data_4["Normal duct","WT"]
data_df <- data.frame(cell_type = rep(c("Normal cell", "Normal duct cell", "PanIN cell", "PDAC cell"), each = 2),
                      probe_count = c(normal_var_count, normal_ref_count, normal_duct_var_count, normal_duct_ref_count, panin_var_count, panin_ref_count, pdac_var_count, pdac_ref_count),
                      probe_total = rep(c(sum(normal_var_count, normal_ref_count), sum(normal_duct_var_count, normal_duct_ref_count), sum(panin_var_count, panin_ref_count), sum(pdac_var_count, pdac_ref_count)), each = 2),
                      cell_count = rep(paste0(c(n_cell_normal, n_cell_normal_duct, n_cell_panin, n_cell_pdac)," cells"), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 4))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "Normal duct cell","PanIN cell","PDAC cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants) # Normal v PDAC 2.2e-16 # 1.004522e-15
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants) # PanIN v PDAC 0.2521 # 1
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants) # Normal v PanIN 2.773e-15 # 1.164688e-13
prop_test_adj_p_c1_nd <- p.adjust(variant_prop_test_4$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PanIN 1.771e-10 # 7.440172e-09
prop_test_adj_p_c2_nd <- p.adjust(variant_prop_test_5$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PDAC 3.027e-13 # 1.271176e-11
df_max_round = max(data_df$probe_total)
dir.create(paste0(out_dir,"/",sample))
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.6*df_max_round)) +
    # scale_y_continuous(limits = c(1,10000), trans='log10') +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = df_max_round*1.1, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("PDAC cell", "PanIN cell")),
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = df_max_round*1.5, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = df_max_round*1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_nd, format = "E", digits = 2)),
                y_position = df_max_round*1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_nd, format = "E", digits = 2)),
                y_position = df_max_round*1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = probe_total),
              color = "black", nudge_y = max(data_df$probe_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Count of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_count_p_value_alt_PanIN_vs_PDAC_figure5_supplement.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()
normal_total = normal_var_count + normal_ref_count
normal_duct_total = normal_duct_var_count + normal_duct_ref_count
pdac_total = pdac_var_count + pdac_ref_count
panin_total = panin_var_count + panin_ref_count
data_df <- data.frame(cell_type = rep(c("Normal cell", "Normal duct cell", "PanIN cell", "PDAC cell"), each = 2),
                      proportion = c(normal_var_count/normal_total, 
                                     normal_ref_count/normal_total, 
                                     normal_duct_var_count/normal_duct_total, 
                                     normal_duct_ref_count/normal_duct_total,
                                     pdac_var_count/pdac_total, 
                                     pdac_ref_count/pdac_total, 
                                     panin_var_count/panin_total, 
                                     panin_ref_count/panin_total),
                      prop_total = rep(c(1, 1, 1, 1), each = 2),
                      cell_count = rep(paste0(c(n_cell_normal, n_cell_normal_duct, n_cell_panin, n_cell_pdac)," cells"), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 4))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "Normal duct cell", "PanIN cell", "PDAC cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants) # Normal v PDAC 2.2e-16 # 1.004522e-15
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants) # PanIN v PDAC 0.2521 # 1
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants) # Normal v PanIN 2.773e-15 # 1.164688e-13
prop_test_adj_p_c1_nd <- p.adjust(variant_prop_test_4$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PanIN 1.771e-10 # 7.440172e-09
prop_test_adj_p_c2_nd <- p.adjust(variant_prop_test_5$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PDAC 3.027e-13 # 1.271176e-11
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.6), breaks = seq(from = 0, to = 1, by = 0.2)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.1, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("PDAC cell", "PanIN cell")),
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.5, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_nd, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_nd, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_proportion_p_value_alt_PanIN_vs_PDAC_figure5_supplement.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

# Figure 5 supplement and figure 5c
# conda activate seurat5
# plotting the zoomed in images for figure 5 related supplement
# Plotting features for the pre-cancer progression part zoomed in on regions of normal, PanIN, and PDAC
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
i = 26 # HT270P1-S1H1A1US2_1
# showing G12D and G12V at the same time.
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# since the cell types were generated on an object that was subset there are some 
# cells that were removed when we filtered certain genes out of the matrix prior to integration.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
unknown_low_quality_labels_v3 <- unique(unknown_low_quality_labels_v3,"LowCount")
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive")))] <- "Normal duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("Normal","Normal duct","PanIN","Neoplastic","LowCount/Unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-21-A"
G12V_key <- "KRAS-p-G12V-ALT-21-T"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS p.G12 probe not detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS p.G12 reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS p.G12V variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS p.G12D variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS p.G12D p.G12V variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS p.G12 probe not detected","KRAS p.G12 reference","KRAS p.G12V variant","KRAS p.G12D variant","KRAS p.G12D p.G12V variant"))
zoom1.1 <- Crop(obj[["fov"]], x = c(4900, 5300), y = c(2370,2770), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9','#eee461','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9','#eee461',"#CC79A7",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

zoom2.1 <- Crop(obj[["fov"]], x = c(4000, 4400), y = c(1550,1950), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom2.1"]] <- zoom2.1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom2.1"]]) <- "segmentation"
# The colors are based on the order in which the factor occurs in the metadata. 
# Additionally some factor levels are not present in the zoom1.2 fov. 
# I use a differently ordered color vector to ensure that colors are consistent across plots. 
# list("normal" = '#009E73',"normal_duct" = "#56B4E9", "PDAC" = "#CC79A7", "PanIN" = "#eee461", "low_quality/unknown" = "#5c5c5c")
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9','#5c5c5c',"#eee461",'#CC79A7'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2.1_darkblue_ceiling1_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in supplementary figure
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9','#5c5c5c',"#eee461",'#CC79A7'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom2.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom2.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom2.1_darkblue_ceiling1_legend_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in supplementary figure
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

zoom3.1 <- Crop(obj[["fov"]], x = c(1800, 2200), y = c(845,1245), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom3.1"]] <- zoom3.1 # this resolution (1000 µm by 1000 µm) ended up being too low when the plots were shrunk down. 
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
DefaultBoundary(obj[["zoom3.1"]]) <- "segmentation"
# The colors are based on the order in which the factor occurs in the metadata. 
# Additionally some factor levels are not present in the zoom1.2 fov. 
# I use a differently ordered color vector to ensure that colors are consistent across plots. 
# list("normal" = '#009E73',"normal_duct" = "#56B4E9", "PDAC" = "#CC79A7", "PanIN" = "#eee461", "low_quality/unknown" = "#5c5c5c")
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9',"#CC79A7",'#5c5c5c','#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom3.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom3.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom3.1_darkblue_ceiling1_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in supplementary figure
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9',"#CC79A7",'#5c5c5c','#eee461'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="KRAS_detected", cols = c('darkblue','yellow','red',"green","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom3.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom3.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom3.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom3.1_darkblue_ceiling1_legend_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in supplementary figure
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
# The colors are based on the order in which the factor occurs in the metadata. 
# Additionally some factor levels are not present in the zoom1.2 fov. 
# I use a differently ordered color vector to ensure that colors are consistent across plots. 
# list("normal" = '#009E73',"normal_duct" = "#56B4E9", "PDAC" = "#CC79A7", "PanIN" = "#eee461", "low_quality/unknown" = "#5c5c5c")
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9','#eee461',"#CC79A7",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="KRAS_detected", cols = c('darkblue','yellow','green',"red","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in supplementary figure
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="neoplasm_duct_normal_unknown", cols = c('#009E73','#56B4E9','#eee461',"#CC79A7",'#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="KRAS_detected", cols = c('darkblue','yellow','green',"red","#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "fov.with.snvs", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "fov.with.snvs", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_fov.with.snvs_darkblue_ceiling1_legend_fig5_supplement.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) # use this in supplementary figure
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# this is for the stacked bar plot for figure 5 supplement
library(epitools)
library(effectsize)
variant_key <- "KRAS-p-G12D-ALT-21-A"
wt_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
obj$tumor_cell_call <- obj$neoplasm_duct_normal_unknown 
# obj$tumor_cell_call <- factor(obj$tumor_cell_call, levels = c("Normal","PanIN","Cancer_clone_2","low_quality/unknown"))
cell_freq = as.data.frame(table(obj$tumor_cell_call, obj@meta.data[paste0(variant_key,"_",assay,"_call")][[1]]))
n_cell_panin = sum(cell_freq[(cell_freq$Var1 == "PanIN"),"Freq"])
n_cell_pdac = sum(cell_freq[(cell_freq$Var1 == "Neoplastic"),"Freq"])
n_cell_normal = sum(cell_freq[(cell_freq$Var1 == "Normal"),"Freq"])
n_cell_normal_duct = sum(cell_freq[(cell_freq$Var1 == "Normal duct"),"Freq"])
data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PDAC", "Normal"))
data_1["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"), paste0(wt_key,"_",assay,"_count")])
data_1["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"), paste0(variant_key,"_",assay,"_count")])
data_1["PDAC","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Neoplastic"), paste0(wt_key,"_",assay,"_count")])
data_1["PDAC","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Neoplastic"), paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_1 <- chisq.test(data_1, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_1)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_1
# X-squared = 70.424, df = NA, p-value = 0.0004998
x <- data_1[,'ALT']
n <- rowSums(data_1)
variant_prop_test_1 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_1)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 70.424, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.4838234 1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.9896907 0.3921569 
oddsratio.out = oddsratio(as.matrix(data_1), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |           95% CI
# -----------------------------
#     148.80     | [19.18, 1154.44]
phi.effectsize <- phi(as.matrix(data_1), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.68       | [0.55, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PDAC", "PanIN"))
data_2["PanIN","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")])
data_2["PanIN","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")])
data_2["PDAC","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Neoplastic"),paste0(wt_key,"_",assay,"_count")])
data_2["PDAC","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Neoplastic"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_2 <- chisq.test(data_2, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_2)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_2
# X-squared = 0.446, df = NA, p-value = 0.6157
x <- data_2[,'ALT']
n <- rowSums(data_2)
variant_prop_test_2 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_2)
#[1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 0.446, df = 1, p-value = 0.2521
# alternative hypothesis: greater
# 95 percent confidence interval:
#     -0.01868585  1.00000000
# sample estimates:
#     prop 1    prop 2 
# 0.9896907 0.9772727
oddsratio.out = oddsratio(as.matrix(data_2), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |        95% CI
# --------------------------
#     2.23       | [0.20, 25.06]
phi.effectsize <- phi(as.matrix(data_2), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.00       | [0.00, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PanIN", "Normal"))
data_3["Normal","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(wt_key,"_",assay,"_count")])
data_3["Normal","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"),paste0(variant_key,"_",assay,"_count")])
data_3["PanIN","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")])
data_3["PanIN","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_3 <- chisq.test(data_3, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_3)
# # [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_3
# X-squared = 61.056, df = NA, p-value = 0.0004998
x <- data_3[,'ALT']
n <- rowSums(data_3)
variant_prop_test_3 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_3)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 61.056, df = 1, p-value = 2.773e-15
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.4696674 1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.9772727 0.3921569 
oddsratio.out = oddsratio(as.matrix(data_3), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |          95% CI
# ----------------------------
#     66.65      | [14.72, 301.84]
phi.effectsize <- phi(as.matrix(data_3), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.66       | [0.52, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_4 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PanIN", "Normal duct"))
data_4["Normal duct","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal duct"),paste0(wt_key,"_",assay,"_count")])
data_4["Normal duct","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal duct"),paste0(variant_key,"_",assay,"_count")])
data_4["PanIN","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(wt_key,"_",assay,"_count")])
data_4["PanIN","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "PanIN"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_4 <- chisq.test(data_4, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_4)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_4
# X-squared = 39.35, df = NA, p-value = 0.0004998
x <- data_4[,'ALT']
n <- rowSums(data_4)
variant_prop_test_4 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_4)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 39.35, df = 1, p-value = 1.771e-10
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.319524 1.000000
# sample estimates:
#     prop 1    prop 2 
# 0.9772727 0.3750000
oddsratio.out = oddsratio(as.matrix(data_4), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |         95% CI
# ---------------------------
#     71.67      | [9.66, 531.43]
phi.effectsize <- phi(as.matrix(data_4), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.63       | [0.46, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
data_5 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("PDAC", "Normal duct"))
data_5["Normal duct","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal duct"),paste0(wt_key,"_",assay,"_count")])
data_5["Normal duct","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal duct"),paste0(variant_key,"_",assay,"_count")])
data_5["PDAC","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Neoplastic"),paste0(wt_key,"_",assay,"_count")])
data_5["PDAC","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Neoplastic"),paste0(variant_key,"_",assay,"_count")])
variant_chisq_test_5 <- chisq.test(data_5, simulate.p.value=TRUE, correct=FALSE)
print("simulate.p.value=TRUE, correct=FALSE")
print(variant_chisq_test_5)
# [1] "simulate.p.value=TRUE, correct=FALSE"
# 
# Pearson's Chi-squared test with simulated p-value (based on 2000
#         replicates)
# 
# data:  data_5
# X-squared = 51.83, df = NA, p-value = 0.0004998
x <- data_5[,'ALT']
n <- rowSums(data_5)
variant_prop_test_5 <- prop.test(x, n , alternative = "greater", correct = FALSE) 
print("prop_test alternative=greater")
print(variant_prop_test_5)
# [1] "prop_test alternative=greater"
# 
# 2-sample test for equality of proportions without continuity correction
# 
# data:  x out of n
# X-squared = 51.83, df = 1, p-value = 3.027e-13
# alternative hypothesis: greater
# 95 percent confidence interval:
#     0.3326472 1.0000000
# sample estimates:
#     prop 1    prop 2 
# 0.9896907 0.3750000
oddsratio.out = oddsratio(as.matrix(data_5), method = "fisher")
print("Odds ratio")
print(oddsratio.out)
# [1] "Odds ratio"
# Odds ratio |           95% CI
# -----------------------------
#     160.00     | [14.01, 1826.92]
phi.effectsize <- phi(as.matrix(data_5), digits = 3)
print("Effect size (Pearson's phi)")
print(phi.effectsize)
# [1] "Effect size (Pearson's phi)"
# Phi (adj.) |       95% CI
# -------------------------
#     0.70       | [0.53, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].
normal_var_count <- data_1["Normal","ALT"]
normal_ref_count <- data_1["Normal","WT"]
pdac_var_count <- data_1["PDAC","ALT"]
pdac_ref_count <- data_1["PDAC","WT"]
panin_var_count <- data_3["PanIN","ALT"]
panin_ref_count <- data_3["PanIN","WT"]
normal_duct_var_count <- data_4["Normal duct","ALT"]
normal_duct_ref_count <- data_4["Normal duct","WT"]
data_df <- data.frame(cell_type = rep(c("Normal cell", "Normal duct cell", "PanIN cell", "PDAC cell"), each = 2),
                      probe_count = c(normal_var_count, normal_ref_count, normal_duct_var_count, normal_duct_ref_count, panin_var_count, panin_ref_count, pdac_var_count, pdac_ref_count),
                      probe_total = rep(c(sum(normal_var_count, normal_ref_count), sum(normal_duct_var_count, normal_duct_ref_count), sum(panin_var_count, panin_ref_count), sum(pdac_var_count, pdac_ref_count)), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal, n_cell_normal_duct, n_cell_panin, n_cell_pdac)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 4))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "Normal duct cell","PanIN cell","PDAC cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants) # Normal v PDAC 2.2e-16 # 1.004522e-15
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants) # PanIN v PDAC 0.2521 # 1
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants) # Normal v PanIN 2.773e-15 # 1.164688e-13
prop_test_adj_p_c1_nd <- p.adjust(variant_prop_test_4$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PanIN 1.771e-10 # 7.440172e-09
prop_test_adj_p_c2_nd <- p.adjust(variant_prop_test_5$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PDAC 3.027e-13 # 1.271176e-11
df_max_round = max(data_df$probe_total)
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5*df_max_round)) +
    # scale_y_continuous(limits = c(1,10000), trans='log10') +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = df_max_round*1.1, tip_length = 0, vjust=0) +
    # geom_signif(comparisons=list(c("PDAC cell", "PanIN cell")), 
    #             annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
    #             y_position = df_max_round*1.5, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = df_max_round*1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_nd, format = "E", digits = 2)),
                y_position = df_max_round*1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_nd, format = "E", digits = 2)),
                y_position = df_max_round*1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = probe_total),
              color = "black", nudge_y = max(data_df$probe_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Count of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_count_p_value_alt_PanIN_vs_PDAC_figure5_supplement.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()
normal_total = normal_var_count + normal_ref_count
normal_duct_total = normal_duct_var_count + normal_duct_ref_count
pdac_total = pdac_var_count + pdac_ref_count
panin_total = panin_var_count + panin_ref_count
data_df <- data.frame(cell_type = rep(c("Normal cell", "Normal duct cell", "PanIN cell", "PDAC cell"), each = 2),
                      proportion = c(normal_var_count/normal_total, 
                                     normal_ref_count/normal_total, 
                                     normal_duct_var_count/normal_duct_total, 
                                     normal_duct_ref_count/normal_duct_total,
                                     pdac_var_count/pdac_total, 
                                     pdac_ref_count/pdac_total, 
                                     panin_var_count/panin_total, 
                                     panin_ref_count/panin_total),
                      prop_total = rep(c(1, 1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal, n_cell_normal_duct, n_cell_panin, n_cell_pdac)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 4))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "Normal duct cell", "PanIN cell", "PDAC cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants) # Normal v PDAC 2.2e-16 # 1.004522e-15
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants) # PanIN v PDAC 0.2521 # 1
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants) # Normal v PanIN 2.773e-15 # 1.164688e-13
prop_test_adj_p_c1_nd <- p.adjust(variant_prop_test_4$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PanIN 1.771e-10 # 7.440172e-09
prop_test_adj_p_c2_nd <- p.adjust(variant_prop_test_5$p.value, method = "bonferroni", n = number_of_variants) # Normal duct v PDAC 3.027e-13 # 1.271176e-11
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.1, tip_length = 0, vjust=0) +
    # geom_signif(comparisons=list(c("PDAC cell", "PanIN cell")), 
    #             annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
    #             y_position = 1.5, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PanIN cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_nd, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal duct cell", "PDAC cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_nd, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_proportion_p_value_alt_PanIN_vs_PDAC_figure5_supplement.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

# This is for figure 5c
obj$cell_type_reduced <- obj$cell_type
obj$cell_type_reduced[obj$cell_type %in% c("Islet_delta","Islet_INS+GCG+","Islet_gamma","Islet_alpha","Islet_beta")] <- "Islet"
obj$cell_type_reduced[obj$cell_type %in% c("Dendritic_cell","Macrophage")] <- "Myeloid"
obj$cell_type_reduced[obj$cell_type %in% c("CD8+T_cell","CD4+T_cell","Naive_T_cell","NK_cell","Treg")] <- "T/NK cell"

cancer_genes = c("KRT17","KRT7","LAMC2","KRT19") #,"TFF1") excluding TFF1 because it is high in both pre-cancer and cancer and the goal here is to distinguish between them and also normal.
pre_cancer_genes = c("MUC5AC","TFF2","PGC")
normal_duct_genes = c("ONECUT2","CRP","CFTR","SOX9","FXYD2")
# remaking the expression and genotype heatmap plots
obj <- AddModuleScore(obj, features = list(cancer_genes), assay="SCT", name="PDAC_score", ctrl = 15, seed=1234) # you need to give AddModuleScore a list (e.g. list()) instead of a vector (e.g. c()) in order to get a single score for the entire set of genes. Otherwise, the function will generate one score for each gene.
obj <- AddModuleScore(obj, features = list(pre_cancer_genes), assay="SCT", name="PanIN_score", ctrl = 15, seed=1234)
obj <- AddModuleScore(obj, features = list(normal_duct_genes), assay="SCT", name="Normal_duct_score", ctrl = 15, seed=1234)
test <- AggregateExpression(obj, assays = "Xenium.with.snvs", features = c("KRAS-p-G12D-ALT-21-A","KRAS-p-G12V-WT"),group.by="cell_type_reduced")
print(as.data.frame(test$Xenium.with.snvs))
#                      Acinar Adipocyte Duct-like-1 Duct-like-2 Endothelial iCAF
# KRAS-p-G12D-ALT-21-A      5         0           2           1           2    1
# KRAS-p-G12V-WT            6         0           5           0           1    6
#                      Islet Mast myCAF Myeloid PanIN Plasma smooth-muscle
# KRAS-p-G12D-ALT-21-A     0    0     8       2    87      0             0
# KRAS-p-G12V-WT           7    0     0       8     2      0             0
#                      T-NK-cell Tumor Vascular-smooth-muscle
# KRAS-p-G12D-ALT-21-A         1    96                      1
# KRAS-p-G12V-WT               1     1                      2
group_ident = "cell_type_reduced"
aggregate_snvs_result <- AggregateExpression(obj, assays = "Xenium.with.snvs", features = c("KRAS-p-G12D-ALT-21-A","KRAS-p-G12V-WT"), group.by=group_ident)
aggregate_snvs = as.data.frame(aggregate_snvs_result$Xenium.with.snvs)
aggregate_kras_result = AggregateExpression(obj, assays = "Xenium.with.snvs", features = c("KRAS"), group.by=group_ident)
aggregate_kras = as.data.frame(aggregate_kras_result$Xenium.with.snvs)
rownames(aggregate_kras) <- c("KRAS")
colnames(aggregate_snvs) = gsub("-","_",colnames(aggregate_snvs))
colnames(aggregate_kras) = gsub("-","_",colnames(aggregate_kras))
module_avg_df = data.frame(PDAC_score1 = rep(NA,length(colnames(aggregate_snvs))),
                           PanIN_score1 = rep(NA,length(colnames(aggregate_snvs))),
                           Normal_duct_score1 = rep(NA,length(colnames(aggregate_snvs))),
                           row.names = colnames(aggregate_snvs))
module_avg_df = t(module_avg_df)
cell_type_vec = colnames(module_avg_df)
score_vec = rownames(module_avg_df)
for (cell_type in cell_type_vec) {
    for (score in score_vec) {
        cell_scores = obj@meta.data[(obj$cell_type_reduced == cell_type),score]
        cell_score = mean(cell_scores)
        module_avg_df[score,cell_type] = cell_score
    }
} 
aggregate_snvs = rbind(aggregate_kras, aggregate_snvs, module_avg_df)
aggregate_snvs = as.data.frame(t(aggregate_snvs))
for (col_name in colnames(aggregate_snvs)) {
    aggregate_snvs[,col_name] = as.numeric(aggregate_snvs[,col_name])
}
print(aggregate_snvs)
# row_names <- rownames(aggregate_snvs)
# row_names[grepl("Vascular_smooth_muscle",row_names)] <- "vSMC"
# row_names[grepl("smooth_muscle",row_names)] <- "SMC"
# rownames(aggregate_snvs) <- row_names
aggregate_snvs[,"cell_types"] = rownames(aggregate_snvs)
aggregate_snvs$neoplastic_labels = NA
aggregate_snvs$neoplastic_labels[aggregate_snvs$cell_types %in% c("Tumor")] = "PDAC"
aggregate_snvs$neoplastic_labels[aggregate_snvs$cell_types %in% c("PanIN")] = "PanIN"
aggregate_snvs$neoplastic_labels[aggregate_snvs$cell_types %in% c("Ductal_reactive","Duct_like_1","Duct_like_2")] = "Normal duct"
aggregate_snvs$neoplastic_labels[(!(aggregate_snvs$cell_types %in% c("Tumor","PanIN","Ductal_reactive","Duct_like_1","Duct_like_2")))] = "Normal"
aggregate_snvs$neoplastic_labels[aggregate_snvs$cell_types %in% c("LowCount")] = "LowCount/Unknown"
#re-scaling the Aggregate expression to be from 0-1
aggregate_snvs[,"KRAS-p-G12D-ALT-21-A"] <- aggregate_snvs[,"KRAS-p-G12D-ALT-21-A"]/max(aggregate_snvs[,"KRAS-p-G12D-ALT-21-A"])
aggregate_snvs[,"KRAS-p-G12V-WT"] <- aggregate_snvs[,"KRAS-p-G12V-WT"]/max(aggregate_snvs[,"KRAS-p-G12V-WT"])
aggregate_snvs[,"KRAS"] <- aggregate_snvs[,"KRAS"]/max(aggregate_snvs[,"KRAS"])
library(ComplexHeatmap)
library(circlize)
colnames(aggregate_snvs) <- c("KRAS gene","KRAS p.G12D Variant","KRAS p.G12D Reference","PDAC score","PanIN score","Normal duct score", "Cell_types","Simplified_cell_type")
scores_df = aggregate_snvs[,c("Normal duct score","PanIN score","PDAC score")]
KRAS_df = aggregate_snvs[,c("KRAS p.G12D Reference","KRAS p.G12D Variant")]
kras_gene_df = aggregate_snvs[,c("KRAS gene")]
neoplastic_label_anno = factor(aggregate_snvs[,c("Simplified_cell_type")], levels = c("Normal","Normal duct","PanIN","PDAC","LowCount/Unknown"))
scores_df = t(scores_df)
KRAS_df = t(KRAS_df)
kras_gene_df = t(kras_gene_df)
colnames(kras_gene_df) = rownames(aggregate_snvs)
rownames(kras_gene_df) = colnames(aggregate_snvs)[grepl("gene", colnames(aggregate_snvs))]
# hcl color palettes hcl_palettes for heatmaps: https://colorspace.r-forge.r-project.org/reference/hcl_palettes.html
column_anno = HeatmapAnnotation(Simplified_cell_type = neoplastic_label_anno, 
                                col = list(Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7", "LowCount/Unknown" = "#5c5c5c")), 
                                na_col = "black",
                                annotation_name_side = "right")
scores_colors = colorRamp2(c(min(scores_df), mean(c(min(scores_df),max(scores_df))), max(scores_df)), c("blue", "white", "red"))
KRAS_colors = colorRamp2(c(0, 0.5, 1), hcl_palette = "Viridis") #"#4B0055FF" "#009B95FF" "#FDE333FF"
kras_gene_colors = colorRamp2(c(0, 0.5, 1), hcl_palette = "Inferno") #"#040404FF" "#C53270FF" "#FFFE9EFF"
pdf(paste0(out_dir,"/",sample,"/",sample,"_KRAS_scores_heatmap.pdf"),useDingbats = F, width=10,height=3.3)
scores = Heatmap(scores_df, 
                 name = "Cell type score", 
                 col = scores_colors,
                 show_column_dend = FALSE,
                 show_row_dend = FALSE,
                 row_order = c("Normal duct score","PanIN score","PDAC score"),
                 column_order = c('LowCount','Mast','Endothelial','Myeloid','T/NK cell','Adipocyte','Glial','SMC',"B_cell",'Plasma','iCAF','Acinar','Islet','myCAF','Duct_like_1','Duct_like_2',"Ductal_reactive",'PanIN','Tumor'),
                 row_names_side = "right")
KRAS = Heatmap(KRAS_df,
               name="Scaled aggregate\nallele expression",
               col = KRAS_colors,
               show_column_dend = FALSE,
               show_row_dend = FALSE,
               row_order = c("KRAS p.G12D Variant","KRAS p.G12D Reference"),
               column_order = c('LowCount','Mast','Endothelial','Myeloid','T/NK cell','Adipocyte','Glial','SMC',"B_cell",'Plasma','iCAF','Acinar','Islet','myCAF','Duct_like_1','Duct_like_2',"Ductal_reactive",'PanIN','Tumor'),
               row_names_side = "right")
KRAS_gene = Heatmap(kras_gene_df,
                    name="Scaled aggregate\nraw gene expression",
                    col = kras_gene_colors,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    row_order = c("KRAS gene"),
                    column_order = c('LowCount','Mast','Endothelial','Myeloid','T/NK cell','Adipocyte','Glial','SMC',"B_cell",'Plasma','iCAF','Acinar','Islet','myCAF','Duct_like_1','Duct_like_2',"Ductal_reactive",'PanIN','Tumor'),
                    row_names_side = "right",
                    top_annotation = column_anno)
draw(KRAS_gene)
draw(KRAS)
draw(scores)

combined = KRAS_gene %v% KRAS %v% scores
draw(combined) # this is figure 5c
dev.off()
#making the heatmap from genes as well as module scores.
cancer_genes = c("KRT17","KRT7","LAMC2","KRT19") #,"TFF1") excluding TFF1 because it is high in both pre-cancer and cancer and the goal here is to distinguish between them and also normal.
pre_cancer_genes = c("MUC5AC","TFF2","PGC")
normal_duct_genes = c("ONECUT2","CRP","CFTR","SOX9","FXYD2")
aggregate_snvs_result <- AggregateExpression(obj, assays = "Xenium.with.snvs", features = c("KRAS-p-G12D-ALT-21-A","KRAS-p-G12V-WT"), group.by=group_ident)
aggregate_genes_result <- AggregateExpression(obj, assays = "SCT", features = c(normal_duct_genes, pre_cancer_genes, cancer_genes), group.by=group_ident)
aggregate_snvs_v2 = as.data.frame(aggregate_snvs_result$Xenium.with.snvs)
aggregate_genes_v2 = as.data.frame(aggregate_genes_result$SCT)
v2_snvs_genes = rbind(aggregate_snvs_v2, log10(aggregate_genes_v2))
v2_snvs_genes = as.data.frame(t(v2_snvs_genes))
for (col_name in colnames(v2_snvs_genes)) {
    v2_snvs_genes[,col_name] = as.numeric(v2_snvs_genes[,col_name])
}
v2_snvs_genes[,"cell_types"] = rownames(v2_snvs_genes)
v2_snvs_genes$cell_types = gsub("-","_",v2_snvs_genes$cell_types)
# row_names <- rownames(v2_snvs_genes)
rownames(v2_snvs_genes) <- v2_snvs_genes$cell_types
v2_snvs_genes$neoplastic_labels = NA
v2_snvs_genes$neoplastic_labels[v2_snvs_genes$cell_types %in% c("Tumor")] = "PDAC"
v2_snvs_genes$neoplastic_labels[v2_snvs_genes$cell_types %in% c("PanIN")] = "PanIN"
v2_snvs_genes$neoplastic_labels[v2_snvs_genes$cell_types %in% c("Duct_like_1","Duct_like_2")] = "Normal duct"
v2_snvs_genes$neoplastic_labels[(!(v2_snvs_genes$cell_types %in% c("Tumor","PanIN","Duct_like_1","Duct_like_2")))] = "Normal"
v2_snvs_genes$neoplastic_labels[v2_snvs_genes$cell_types %in% c("LowCount")] = "LowCount/Unknown"
v2_snvs_genes[,"KRAS-p-G12D-ALT-21-A"] <- v2_snvs_genes[,"KRAS-p-G12D-ALT-21-A"]/max(v2_snvs_genes[,"KRAS-p-G12D-ALT-21-A"])
v2_snvs_genes[,"KRAS-p-G12V-WT"] <- v2_snvs_genes[,"KRAS-p-G12V-WT"]/max(v2_snvs_genes[,"KRAS-p-G12V-WT"])
col_names <- colnames(v2_snvs_genes)
col_names[grepl("KRAS-p-G12D-ALT-21-A", col_names)] <- "KRAS p.G12D Variant"
col_names[grepl("KRAS-p-G12V-WT", col_names)] <- "KRAS p.G12D Reference"
col_names[grepl("cell_types", col_names)] <- "Cell_types"
col_names[grepl("neoplastic_labels", col_names)] <- "Simplified_cell_type"
colnames(v2_snvs_genes) <- col_names
for (colname in colnames(v2_snvs_genes)) {
    v2_snvs_genes[,colname][is.infinite(v2_snvs_genes[,colname])] <- NA
}
scores_df = v2_snvs_genes[,(!(grepl("KRAS",colnames(v2_snvs_genes))) & !(grepl("Cell",colnames(v2_snvs_genes))) & !(grepl("Simplified",colnames(v2_snvs_genes))))]
KRAS_df = v2_snvs_genes[,grepl("KRAS",colnames(v2_snvs_genes))]
neoplastic_label_anno = factor(v2_snvs_genes[,c("Simplified_cell_type")], levels = c("Normal","Normal duct","PanIN","PDAC","LowCount/Unknown"))
scores_df = t(scores_df)
KRAS_df = t(KRAS_df)
scores_df[is.infinite(scores_df)] = NA # infinite values are the result of log10(0). By setting them to NA we can set their color to the lowest color on the color break scale which is 0.
KRAS_df[is.infinite(KRAS_df)] = NA # infinite values are the result of log10(0). By setting them to NA we can set their color to the lowest color on the color break scale which is 0.
column_anno = HeatmapAnnotation(Simplified_cell_type = neoplastic_label_anno, 
                                col = list(Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7","LowCount/Unknown" = "#5c5c5c")), 
                                na_col = "black",
                                annotation_name_side = "right")
scores_colors = colorRamp2(c(0, mean(c(0,max(scores_df[!(is.na(scores_df))]))), max(scores_df[!(is.na(scores_df))])), c("blue", "white", "red"))
#KRAS_colors = colorRamp2(c(0, mean(c(0,max(KRAS_df[!(is.na(KRAS_df))]))), max(KRAS_df[!(is.na(KRAS_df))])), c("blue", "white", "red"))
KRAS_colors = colorRamp2(c(0, mean(c(0,max(KRAS_df[!(is.na(KRAS_df))]))), max(KRAS_df[!(is.na(KRAS_df))])), hcl_palette = "Viridis")
# attr(,"breaks")
# [1] 0.0000000 0.9911356 1.9822712
# attr(,"colors")
# [1] "#4B0055FF" "#009B95FF" "#FDE333FF" 
pdf(paste0(out_dir,"/",sample,"/",sample,"_KRAS_genes_heatmap.pdf"),useDingbats = F, width=9,height=5)
scores = Heatmap(scores_df, 
                 name = "log10(Gene_aggregate_Expression)", 
                 col = scores_colors,
                 #show_column_dend = FALSE,
                 #show_row_dend = FALSE,
                 #row_order = c(normal_duct_genes, pre_cancer_genes, cancer_genes),
                 #row_order = c("CRP","FXYD2","CFTR","SOX9","","","","",'','','','')
                 row_order = c("FXYD2","CFTR","CRP","ONECUT2","SOX9","TFF2","PGC","MUC5AC",'KRT19','KRT7','LAMC2',"KRT17"),
                 #column_order = c("Tumor",'PanIN','Duct_like_2',"Duct_like_1","Acinar","myCAF","Islet","Myeloid","iCAF","T_NK_cell","Endothelial","Vascular_smooth_muscle","Mast","smooth_muscle","Plasma","Adipocyte"),
                 #column_order = c('','Endothelial','Myeloid','T_NK_cell','Adipocyte','Vascular_smooth_muscle','smooth_muscle','Plasma','iCAF','Acinar','Islet','myCAF','Duct_like_1','Duct_like_2','PanIN','Tumor'),
                 column_order = c('LowCount','Mast','Endothelial','Myeloid','T/NK cell','Adipocyte','Glial','SMC',"B_cell",'Plasma','iCAF','Acinar','Islet','myCAF','Duct_like_1','Duct_like_2',"Ductal_reactive",'PanIN','Tumor'),
                 row_names_side = "right",
                 row_split = factor(c(rep("Ductal",5), rep("PanIN",3), rep("PDAC",4)), levels = c("Ductal","PanIN","PDAC")),
                 na_col = "blue")
KRAS = Heatmap(KRAS_df,
               name="Scaled aggregate\nallele expression",
               col = KRAS_colors,
               show_column_dend = FALSE,
               show_row_dend = FALSE,
               row_order = c("KRAS p.G12D Variant","KRAS p.G12D Reference"),
               #column_order = c("Tumor",'PanIN','Duct_like_2',"Duct_like_1","Acinar","myCAF","Islet","Myeloid","iCAF","T_NK_cell","Endothelial","Vascular_smooth_muscle","Mast","smooth_muscle","Plasma","Adipocyte"),
               column_order = c('LowCount','Mast','Endothelial','Myeloid','T/NK cell','Adipocyte','Glial','SMC',"B_cell",'Plasma','iCAF','Acinar','Islet','myCAF','Duct_like_1','Duct_like_2',"Ductal_reactive",'PanIN','Tumor'),
               row_names_side = "right",
               top_annotation = column_anno,
               na_col = "#4B0055FF")
draw(KRAS)
draw(scores)
combined = KRAS %v% scores
draw(combined) # this is supplement
dev.off()
# Following the initial morph layer calculation of HT270P1 I spoke with Li about how we should handle the population of PanIN cells that is detached from all other parts of the tissue in the image. 
# She advised that since it was detached and we cannot be certain from where on the tissue it originated we should exclude it from the analysis. As a result I have circled the tumor region and the nearby PanIN structure in the Xenium explorer data file.
# I will then use the cell barcodes there to re-label just the PanIN cells so that they are excluded form the layer calculation. 
#conda activate seurat5
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/HT270P1-S1H1A1US2_1/detached_morph_blacklist_for_panin/
library(tidyverse)
set.seed(1234)
cell_types <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/HT270P1-S1H1A1US2_1/Manual_individual_cell_type_v1.tsv", sep='\t',header=T)
detached_piece_df <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/HT270P1-S1H1A1US2_1/detached_morph_blacklist_for_panin/Selection_1_cells_stats_fix_header.csv",sep=",",header=T)
PanIN_detached_piece = detached_piece_df$Cell.ID[detached_piece_df$Cluster == "PanIN"]
cell_types$Manual_individual_cell_types_v1[cell_types$cell_id %in% PanIN_detached_piece] = "PanIN_detached_blacklist"
write.table(cell_types,file="/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/manual_individual_cell_types_v1/HT270P1-S1H1A1US2_1/detached_morph_blacklist_for_panin/Manual_individual_cell_type_v1_PanIN_detached_blacklist.tsv", quote=F, sep="\t", col.names=T, row.names=F)
# Fixing the colors of the binned stacked bar plot from Andre so that it is consistent across both plots after he removed the detached piece
# the tab20 colormap from matplotlib (pulled from here: https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2F09509942-4112-474f-9dfb-2698e6a1e4c2%2F26e4e036-74d6-47ba-85f0-6bbd91e51f58%2Ffiles%2Ffunctions%2Fmatplotlib%2Ftab20.m&embed=web)
# it also matches the hex values I get from adobe illustrator when I use the color picker tool
library(reshape2)
panin_barcode_layer_df <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT270P1-S1H1A1US2_1/v7/HT270P1-S1H1A1US2_1_cell_layer_panin_source.csv",sep=",",header=T)
colnames(panin_barcode_layer_df) <- c("barcode","layer")
panin_barcode_layer_df$layer_floor3 <- panin_barcode_layer_df$layer %/% 3 # use floor division to group every 3 layers together for plotting purposes.
panin_barcode_layer_df$layer_bin <- paste0((panin_barcode_layer_df$layer_floor3*30),"-",((panin_barcode_layer_df$layer_floor3+1)*30)) #generate the bins for plotting
rownames(panin_barcode_layer_df) <- panin_barcode_layer_df$barcode
panin_barcode_layer_df$barcode <- NULL
obj <- AddMetaData(obj, panin_barcode_layer_df)
#Layer 0 means that the cell is in one of the tiles that covers the Tumor/panin starting microregion. So for all cells with layer = 0 we will set the layer_bin to NA
obj$layer_bin[obj$layer == 0] <- NA
# when we blacklisted the detached tissue from the morph analysis this was only from the microregion identification step
# the black listed cells were still present in the layer calculation so we comment out the PanIN cells here.
panin_count <- as.data.frame(table(obj$cell_type_reduced, obj$layer_bin))
colnames(panin_count) <- c("cell_type","layer","count")
panin_count$count[(panin_count$cell_type == "PanIN")] <- 0
layer_vec <- unique(panin_count$layer)
panin_count$proportion <- NA
for (layer in layer_vec) {
    total_counts <- sum(panin_count$count[panin_count$layer == layer])
    panin_count$proportion[panin_count$layer == layer] <- panin_count$count[panin_count$layer == layer]/total_counts
}
layer_start_tmp = as.numeric(unlist(strsplit(as.character(panin_count$layer), "-")))
layer_start_vector <- layer_start_tmp[seq(1, length(layer_start_tmp), 2)]
panin_count$layer_start <- layer_start_vector
cell_order <- c("Tumor","PanIN","Ductal_reactive","Duct_like_2",
                "Duct_like_1","Acinar","Islet",
                "Endothelial","SMC","Glial",
                "myCAF","iCAF","Adipocyte",
                "Mast","Plasma","B_cell","Myeloid",
                "T/NK cell","LowCount")
panin_count$sort <- 1:length(panin_count$cell_type)
sorting_vector <- c()
for (layer in sort(unique(panin_count$layer_start))) {
    for (cell in cell_order) {
        sorting_vector <- c(sorting_vector, panin_count$sort[(panin_count$cell_type == cell) & (panin_count$layer_start == layer)])
    }
}
sorted_panin_count <- panin_count[match(sorting_vector, panin_count$sort),]
rownames(sorted_panin_count) <- 1:length(sorted_panin_count$sort)
# rep_cell_order <- rep(cell_order, each = dim(panin_count)[1]/length(cell_order))
# test <- panin_count[match(rep_cell_order, panin_count$cell_type),]
sorted_panin_count$cell_type <- factor(sorted_panin_count$cell_type, levels = cell_order)
tab20_colors <- c('#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2','#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5')
tab20_colors_reorder <- c('#e377c2','#bcbd22','#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#17becf','#9edae5','#f7b6d2','#dbdb8d','#7f7f7f','#c7c7c7')
layer_labels <- unique(sorted_panin_count$layer[order(sorted_panin_count$layer_start)])
sorted_panin_count_plot <- sorted_panin_count[(sorted_panin_count$layer_start < 300),]
p1 <- ggplot(sorted_panin_count_plot, aes(x=layer_start, y = proportion, fill = cell_type)) +
    geom_bar(position="fill", stat="identity") + 
    theme_classic() + 
    scale_x_continuous(breaks=seq(1,300,30),labels=layer_labels[1:10]) +
    scale_fill_manual( breaks = cell_order,
                       values = tab20_colors_reorder[1:19]) +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
morph_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT270P1-S1H1A1US2_1/v7"
pdf(paste0(morph_dir,"/",sample,"_Panin_source_panin_cells_detach_blacklist_stacked_bar_plot_proportion_layers.pdf"),useDingbats = F, height = 5, width = 6) # figure 5d left
print(p1)
dev.off()
write.table(sorted_panin_count,paste0(morph_dir,"/",sample,"_Panin_source_panin_cells_detach_blacklist_stacked_bar_plot_proportion_layers.tsv",sep='\t',quote=F))

tumor_barcode_layer_df <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT270P1-S1H1A1US2_1/v7/HT270P1-S1H1A1US2_1_cell_layer_tumor_1_source.csv",sep=",",header=T)
colnames(tumor_barcode_layer_df) <- c("barcode","layer")
tumor_barcode_layer_df$layer_floor3 <- tumor_barcode_layer_df$layer %/% 3 # use floor division to group every 3 layers together for plotting purposes.
tumor_barcode_layer_df$layer_bin <- paste0((tumor_barcode_layer_df$layer_floor3*30),"-",((tumor_barcode_layer_df$layer_floor3+1)*30)) #generate the bins for plotting
rownames(tumor_barcode_layer_df) <- tumor_barcode_layer_df$barcode
tumor_barcode_layer_df$barcode <- NULL
obj <- AddMetaData(obj, tumor_barcode_layer_df)
# Layer 0 means that the cell is in one of the tiles that covers the Tumor/panin starting microregion. So for all cells with layer = 0 we will set the layer_bin to NA
obj$layer_bin[obj$layer == 0] <- NA
# when we blacklisted the detached tissue from the morph analysis this was only from the microregion identification step
# the black listed cells were still present in the layer calculation so we comment out the PanIN cells here.
tumor_count <- as.data.frame(table(obj$cell_type_reduced, obj$layer_bin))
colnames(tumor_count) <- c("cell_type","layer","count")
tumor_count$count[(tumor_count$cell_type == "Tumor")] <- 0
layer_vec <- unique(tumor_count$layer)
tumor_count$proportion <- NA
for (layer in layer_vec) {
    total_counts <- sum(tumor_count$count[tumor_count$layer == layer])
    tumor_count$proportion[tumor_count$layer == layer] <- tumor_count$count[tumor_count$layer == layer]/total_counts
}
layer_start_tmp = as.numeric(unlist(strsplit(as.character(tumor_count$layer), "-")))
layer_start_vector <- layer_start_tmp[seq(1, length(layer_start_tmp), 2)]
tumor_count$layer_start <- layer_start_vector
cell_order <- c("Tumor","PanIN","Ductal_reactive","Duct_like_2",
                "Duct_like_1","Acinar","Islet",
                "Endothelial","SMC","Glial",
                "myCAF","iCAF","Adipocyte",
                "Mast","Plasma","B_cell","Myeloid",
                "T/NK cell","LowCount")
tumor_count$sort <- 1:length(tumor_count$cell_type)
sorting_vector <- c()
for (layer in sort(unique(tumor_count$layer_start))) {
    for (cell in cell_order) {
        sorting_vector <- c(sorting_vector, tumor_count$sort[(tumor_count$cell_type == cell) & (tumor_count$layer_start == layer)])
    }
}
sorted_tumor_count <- tumor_count[match(sorting_vector, tumor_count$sort),]
rownames(sorted_tumor_count) <- 1:length(sorted_tumor_count$sort)
# rep_cell_order <- rep(cell_order, each = dim(tumor_count)[1]/length(cell_order))
# test <- tumor_count[match(rep_cell_order, tumor_count$cell_type),]
sorted_tumor_count$cell_type <- factor(sorted_tumor_count$cell_type, levels = cell_order)
tab20_colors <- c('#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2','#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5')
tab20_colors_reorder <- c('#e377c2','#bcbd22','#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#17becf','#9edae5','#f7b6d2','#dbdb8d','#7f7f7f','#c7c7c7')
layer_labels <- unique(sorted_tumor_count$layer[order(sorted_tumor_count$layer_start)])
sorted_tumor_count_plot <- sorted_tumor_count[(sorted_tumor_count$layer_start < 300),]
p1 <- ggplot(sorted_tumor_count_plot, aes(x=layer_start, y = proportion, fill = cell_type)) +
    geom_bar(position="fill", stat="identity") + 
    theme_classic() + 
    scale_x_continuous(breaks=seq(1,300,30),labels=layer_labels[1:10]) +
    scale_fill_manual( breaks = cell_order,
                       values = tab20_colors_reorder[1:19]) +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
morph_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/morph/binned_stacked_bar/HT270P1-S1H1A1US2_1/v7"
pdf(paste0(morph_dir,"/",sample,"_tumor_source_panin_cells_detach_blacklist_stacked_bar_plot_proportion_layers.pdf"),useDingbats = F, height = 5, width = 6) # figure 5d right
print(p1)
dev.off()
write.table(sorted_tumor_count,paste0(morph_dir,"/",sample,"_tumor_source_panin_cells_detach_blacklist_stacked_bar_plot_proportion_layers.tsv"),sep='\t',quote=F)

# Figure 6 spatial inset plots:
# i = 19 # SP001P1-Fp1U1 # G12D 
# i = 26 # HT270P1-S1H1A1US2_1 # G12D # replace with G12V HT125P1-S1H8A1U1
# i = 27 # SP002C1-Fp1U2 # G12V
# i = 34 # HT179C1-T1Fp3L5U1 # G12D
# HT242P1
# HT227P1
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
options(future.globals.maxSize= +Inf) # need to increase for SP001P1 #+Inf is 100% overkill but I don't care.
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
i = 19 # SP001P1-Fp1U1 # G12D
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
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
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive")))] <- "Normal duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("Normal","Normal duct","PanIN","PDAC","LowCount/Unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS p.G12 probe not detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS p.G12 reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS p.G12V variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS p.G12D variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS p.G12D p.G12V variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS p.G12 probe not detected","KRAS p.G12 reference","KRAS p.G12V variant","KRAS p.G12D variant","KRAS p.G12D p.G12V variant"))
zoom1.1 <- Crop(obj[["fov"]], x = c(5500, 6000), y = c(5300, 5800), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#FF7F00','#F781BF',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#FF7F00','#F781BF',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
# i = 26 # HT270P1-S1H1A1US2_1 # G12D
# replace with G12V
i = 28 # HT060P1-S1R1Fp1U1 # G12V
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# since the cell types were generated on an object that was subset there are some 
# cells that were removed when we filtered certain genes out of the matrix prior to integration.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
unknown_low_quality_labels_v3 <- unique(unknown_low_quality_labels_v3,"LowCount")
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive")))] <- "Normal duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("Normal","Normal duct","PanIN","Neoplastic","LowCount/Unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS p.G12 probe not detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS p.G12 reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS p.G12V variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS p.G12D variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS p.G12D p.G12V variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS p.G12 probe not detected","KRAS p.G12 reference","KRAS p.G12V variant","KRAS p.G12D variant","KRAS p.G12D p.G12V variant"))
zoom1.1 <- Crop(obj[["fov"]], x = c(11150, 11550), y = c(3890, 4290), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#F781BF','#FF7F00',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#F781BF','#FF7F00',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
i = 32 # HT125P1-S1H8A1U1 # G12V
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# since the cell types were generated on an object that was subset there are some 
# cells that were removed when we filtered certain genes out of the matrix prior to integration.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
unknown_low_quality_labels_v3 <- unique(unknown_low_quality_labels_v3,"LowCount")
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive")))] <- "Normal duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("Normal","Normal duct","PanIN","Neoplastic","LowCount/Unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS p.G12 probe not detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS p.G12 reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS p.G12V variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS p.G12D variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS p.G12D p.G12V variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS p.G12 probe not detected","KRAS p.G12 reference","KRAS p.G12V variant","KRAS p.G12D variant","KRAS p.G12D p.G12V variant"))
# use zoom1.1
zoom1.1 <- Crop(obj[["fov"]], x = c(1600, 2000), y = c(3750, 4150), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#F781BF','#FF7F00',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#F781BF','#FF7F00',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
i = 27 # SP002C1-Fp1U2 # G12V
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# since the cell types were generated on an object that was subset there are some 
# cells that were removed when we filtered certain genes out of the matrix prior to integration.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
unknown_low_quality_labels_v3 <- unique(unknown_low_quality_labels_v3,"LowCount")
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive")))] <- "Normal duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("Normal","Normal duct","PanIN","Neoplastic","LowCount/Unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-21-A"
G12V_key <- "KRAS-p-G12V-ALT-21-T"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS p.G12 probe not detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS p.G12 reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS p.G12V variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS p.G12D variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS p.G12D p.G12V variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS p.G12 probe not detected","KRAS p.G12 reference","KRAS p.G12V variant","KRAS p.G12D variant","KRAS p.G12D p.G12V variant"))
# use zoom1.1
zoom1.1 <- Crop(obj[["fov"]], x = c(4650, 5150), y = c(2350, 2850), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#F781BF','#FF7F00',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#F781BF','#FF7F00',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
i = 33 # HT179C1-T1Fp3L5U1 # G12D
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
DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
Idents(obj) <- "seurat_clusters"
DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
obj$barcode <- rownames(obj@meta.data)
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
# since the cell types were generated on an object that was subset there are some 
# cells that were removed when we filtered certain genes out of the matrix prior to integration.
obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
unknown_low_quality_labels_v3 <- unique(unknown_low_quality_labels_v3,"LowCount")
obj$neoplasm_normal_unknown <- NA
obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
#obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
#obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
#obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
obj$neoplasm_duct_normal_unknown <- NA
obj$neoplasm_duct_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
obj$neoplasm_duct_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
obj$neoplasm_duct_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
obj$neoplasm_duct_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2","Ductal_reactive")))] <- "Normal duct"
obj$neoplasm_duct_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
obj$neoplasm_duct_normal_unknown <- factor(obj$neoplasm_duct_normal_unknown, levels = c("Normal","Normal duct","PanIN","Neoplastic","LowCount/Unknown"))
tmp_df <- obj@meta.data[ , c("barcode", "neoplasm_duct_normal_unknown")]
colnames(tmp_df) = c("cell_id","group")
write.table(tmp_df, paste0(out_dir,"/",sample,"/",sample,"_neoplastic_duct_normal_unknown_labels.csv"),sep=',',quote=F)
# there are no normal duct or panin cells in this case
obj$uniform_background_color = "darkblue"
G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
assay = "Xenium.with.snvs"
assay_class = class(obj[[assay]])
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
G12D_tmp <- as.data.frame(counts_df[ ,G12D_key],row.names=rownames(counts_df))
G12V_tmp <- as.data.frame(counts_df[ ,G12V_key],row.names=rownames(counts_df))
ref_tmp <- as.data.frame(counts_df[ ,ref_key], row.names=rownames(counts_df))
obj <- AddMetaData(object = obj, metadata = G12D_tmp, col.name = paste0(G12D_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = G12V_tmp, col.name = paste0(G12V_key,"_",assay,"_count"))
obj <- AddMetaData(object = obj, metadata = ref_tmp, col.name = paste0(ref_key,"_",assay,"_count"))
obj@meta.data["KRAS_detected"] <- "KRAS p.G12 probe not detected"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(ref_key,"_",assay,"_count")] > 0] <- "KRAS p.G12 reference"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12V_key,"_",assay,"_count")] > 0] <- "KRAS p.G12V variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0] <- "KRAS p.G12D variant"
obj@meta.data["KRAS_detected"][obj@meta.data[paste0(G12D_key,"_",assay,"_count")] > 0 & obj@meta.data[paste0(G12V_key,"_",assay,"_count")] ] <- "KRAS p.G12D p.G12V variant"
obj$KRAS_detected <- factor(obj$KRAS_detected, levels = c("KRAS p.G12 probe not detected","KRAS p.G12 reference","KRAS p.G12V variant","KRAS p.G12D variant","KRAS p.G12D p.G12V variant"))
# use zoom1.1
zoom1.1 <- Crop(obj[["fov"]], x = c(5330, 5730), y = c(1550, 1950), coords = c("plot","tissue")) # the x and y coordinates here are the transposed coordinates from the region of interest targeted for the crop because the seurat developers can't keep their x and y coordinates straight for whatever reason.
obj[["zoom1.1"]] <- zoom1.1
#DefaultFOV(obj, assay='Xenium.with.snvs') <- 'zoom1'
#Simplified_cell_type = c("Normal" = "#009E73","Normal duct" = "#56B4E9","PanIN" = "#eee461","PDAC" = "#CC79A7")
DefaultBoundary(obj[["zoom1.1"]]) <- "segmentation"
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#FF7F00','#F781BF',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
# ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#140152", mid = "#BA50DD", high = "#F20089", midpoint = 0.5)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()
p13 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="neoplasm_normal_unknown", cols = c('#009E73','#CC79A7','#5c5c5c'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p14 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="cell_type", border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="KRAS_detected", cols = c('darkblue',"#4DAF4A",'#FF7F00','#F781BF',"#5c5c5c"), border.color=NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
list_of_plots <- list(p13)
list_of_plots[[ length(list_of_plots) + 1 ]] <- p14
list_of_plots[[ length(list_of_plots) + 1 ]] <- p15
jiterables = length(feature_list)
offset = length(list_of_plots)
DefaultAssay(obj) <- "Xenium.with.snvs"
for (j in 1:jiterables) {
    p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = feature_list[j], min.cutoff = 0, max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.color= NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    list_of_plots[[j+offset]] <- p15
}
DefaultAssay(obj) <- "SCT"
p15 <- rasterize(ImageFeaturePlot(obj, fov = "zoom1.1", features = "KRAS", size = 0.15, min.cutoff = 0, max.cutoff = 1, scale = "feature", border.color= NA, axes = TRUE) + scale_fill_continuous(type = "viridis") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title =element_blank()) + coord_fixed(ratio=1) , layers="Polygon", dpi=600) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
p15 <- rasterize(ImageDimPlot(obj, fov = "zoom1.1", group.by="uniform_background_color", cols = c('darkblue'), border.color = NA, axes = TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers="Polygon", dpi=600)
offset = length(list_of_plots)
list_of_plots[[1+offset]] <- p15
pdf(paste0(out_dir,"/",sample,"/",sample,"_neoplastic_normal_unknown_v3_image_2color_zoom1.1_darkblue_ceiling1_legend_fig6.pdf"),useDingbats = F,height=5, width = 5*(length(list_of_plots))) #this might be barely okay but it doesn't really offer much benefit beyond what we get from the Xenium explorer images
print(ggarrange(plotlist = list_of_plots, ncol = length(list_of_plots), nrow = 1))
dev.off()

# plotting the stacked bar plots showing proportion of KRAS p.G12D alleles in normal cells 
library(epitools)
library(effectsize)
G12D_sections <- c('HT227P1-S1H1L1U1',
                   'HT242P1-S1H4L4U1',
                   'SP001P1-Fp1U1',
                   'HT270P1-S1H1A1US2_1',
                   'HT061P1-S1P1A1L1U1',
                   'HT061P1-S1P1A1L4U1',
                   'HT179C1-T1Fp3L5U1')
G12V_sections <- c('SP002C1-Fp1U2',
                   'HT060P1-S1R1Fp1U1',
                   'HT125P1-S1H4A1L1U1',
                   'HT125P1-S1H8A1U1',
                   'HT185P1-S1H2L1U1')
KRAS_reference_sections <- unique(all_sample_summary$sample.ID[!(all_sample_summary$sample.ID %in% G12D_sections) & !(all_sample_summary$sample.ID %in% G12V_sections)])
all_sample_summary_cells <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)

G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"

normal_cell_VAR_G12D_count <- sum(all_sample_summary$Number.of.variant.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)]) +
    sum(all_sample_summary$Number.of.variant.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)]) + 
    sum(all_sample_summary$Number.of.variant.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)]) +
    sum(all_sample_summary$Number.of.variant.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)])
    
normal_cell_REF_G12D_count <- sum(all_sample_summary$Number.of.reference.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)]) + 
                        sum(all_sample_summary$Number.of.reference.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)]) + 
                        sum(all_sample_summary$Number.of.reference.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)]) +
                        sum(all_sample_summary$Number.of.reference.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)])

G12D_cell_VAR_G12D_count <- sum(all_sample_summary$Number.of.variant.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)])
G12D_cell_REF_G12D_count <- sum(all_sample_summary$Number.of.reference.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)])

G12V_cell_VAR_G12D_count <- sum(all_sample_summary$Number.of.variant.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)])
G12V_cell_REF_G12D_count <- sum(all_sample_summary$Number.of.reference.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12D_key)])

normal_cell_VAR_G12D <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)]) +
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)]) + 
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)]) +
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)])

normal_cell_REF_G12D <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)]) + 
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)]) + 
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)]) +
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)])

G12D_cell_VAR_G12D <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)])
G12D_cell_REF_G12D <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)])

G12V_cell_VAR_G12D <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)])
G12V_cell_REF_G12D <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12D_key)])

data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("G12D_cancer_cell", "non-G12D_cancer_and_Normal_cell"))

data_1["non-G12D_cancer_and_Normal_cell","WT"] <- normal_cell_REF_G12D_count
data_1["non-G12D_cancer_and_Normal_cell","ALT"] <- normal_cell_VAR_G12D_count
data_1["G12D_cancer_cell","WT"] <- G12D_cell_REF_G12D_count
data_1["G12D_cancer_cell","ALT"] <- G12D_cell_VAR_G12D_count
x <- data_1[,'ALT']
#print(x)
n <- rowSums(data_1)
#print(n)
variant_prop_test_1 <- prop.test(x, n, alternative = "greater", correct = FALSE)

data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("G12V_cancer_cell", "non-G12D_cancer_and_Normal_cell"))
data_2["non-G12D_cancer_and_Normal_cell","WT"] <- normal_cell_REF_G12D_count
data_2["non-G12D_cancer_and_Normal_cell","ALT"] <- normal_cell_VAR_G12D_count
data_2["G12V_cancer_cell","WT"] <- G12V_cell_REF_G12D_count
data_2["G12V_cancer_cell","ALT"] <- G12V_cell_VAR_G12D_count
x <- data_2[,'ALT']
#print(x)
n <- rowSums(data_2)
#print(n)
variant_prop_test_2 <- prop.test(x, n, alternative = "greater", correct = FALSE) 

data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("G12D_cancer_cell", "G12V_cancer_cell"))
data_3["G12V_cancer_cell","WT"] <- G12V_cell_REF_G12D_count
data_3["G12V_cancer_cell","ALT"] <- G12V_cell_VAR_G12D_count
data_3["G12D_cancer_cell","WT"] <- G12D_cell_REF_G12D_count
data_3["G12D_cancer_cell","ALT"] <- G12D_cell_VAR_G12D_count
x <- data_3[,'ALT']
#print(x)
n <- rowSums(data_3)
#print(n)
variant_prop_test_3 <- prop.test(x, n, alternative = "greater", correct = FALSE) #G12V vs G12D

normal_var_count <- normal_cell_VAR_G12D_count
normal_ref_count <- normal_cell_REF_G12D_count
clone_2_var_count <- G12D_cell_VAR_G12D_count
clone_2_ref_count <- G12D_cell_REF_G12D_count
clone_1_var_count <- G12V_cell_VAR_G12D_count
clone_1_ref_count <- G12V_cell_REF_G12D_count

normal_total = normal_var_count + normal_ref_count
clone_2_total = clone_2_var_count + clone_2_ref_count
clone_1_total = clone_1_var_count + clone_1_ref_count

n_cell_normal = normal_cell_VAR_G12D + normal_cell_REF_G12D
n_cell_clone_2 = G12D_cell_VAR_G12D + G12D_cell_REF_G12D
n_cell_clone_1 = G12V_cell_VAR_G12D + G12V_cell_REF_G12D

data_df <- data.frame(cell_type = rep(c("Normal cell", "G12D cancer cell", "G12V cancer cell"), each = 2),
                      proportion = c(normal_var_count/normal_total,normal_ref_count/normal_total,clone_2_var_count/clone_2_total,clone_2_ref_count/clone_2_total,clone_1_var_count/clone_1_total,clone_1_ref_count/clone_1_total),
                      prop_total = rep(c(1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "G12D cancer cell", "G12V cancer cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "G12D cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("G12D cancer cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of KRAS p.G12D probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/KRAS_p.G12D_barplots_proportion_p_value_alt_g12d_vs_g12v.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

data_df <- data.frame(cell_type = rep(c("Normal cell", "G12D cancer cell", "G12V cancer cell"), each = 2),
                      probe_count = c(normal_var_count,normal_ref_count,clone_2_var_count,clone_2_ref_count,clone_1_var_count,clone_1_ref_count),
                      prop_total = rep(c(n_cell_normal,n_cell_clone_2,n_cell_clone_1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "G12D cancer cell", "G12V cancer cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
df_max_round = max(c(n_cell_normal,n_cell_clone_2,n_cell_clone_1))
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5*df_max_round)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "G12D cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.2*df_max_round, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("G12D cancer cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3*df_max_round, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.4*df_max_round, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of KRAS p.G12D probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/KRAS_p.G12D_barplots_count_p_value_alt_g12d_vs_g12v.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

G12D_key <- "KRAS-p-G12D-ALT-T"
G12V_key <- "KRAS-p-G12V-ALT-A"
ref_key <- "KRAS-p-G12V-WT"
normal_cell_VAR_G12V_count <- sum(all_sample_summary$Number.of.variant.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)]) +
    sum(all_sample_summary$Number.of.variant.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)]) + 
    sum(all_sample_summary$Number.of.variant.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)]) +
    sum(all_sample_summary$Number.of.variant.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)])

normal_cell_REF_G12V_count <- sum(all_sample_summary$Number.of.reference.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)]) + 
    sum(all_sample_summary$Number.of.reference.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% KRAS_reference_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)]) + 
    sum(all_sample_summary$Number.of.reference.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)]) +
    sum(all_sample_summary$Number.of.reference.probes.in.non.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)])

G12D_cell_VAR_G12V_count <- sum(all_sample_summary$Number.of.variant.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)])
G12D_cell_REF_G12V_count <- sum(all_sample_summary$Number.of.reference.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12D_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)])

G12V_cell_VAR_G12V_count <- sum(all_sample_summary$Number.of.variant.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)])
G12V_cell_REF_G12V_count <- sum(all_sample_summary$Number.of.reference.probes.in.cancer.cells[(all_sample_summary$sample.ID %in% G12V_sections) & (all_sample_summary$Variant.probe.name %in% G12V_key)])

normal_cell_VAR_G12V <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)]) +
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)]) + 
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)]) +
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)])

normal_cell_REF_G12V <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)]) + 
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% KRAS_reference_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)]) + 
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)]) +
    sum(all_sample_summary_cells$Number.of.non.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)])

G12D_cell_VAR_G12V <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)])
G12D_cell_REF_G12V <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12D_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)])

G12V_cell_VAR_G12V <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.variant.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)])
G12V_cell_REF_G12V <- sum(all_sample_summary_cells$Number.of.cancer.cells.with.the.reference.probe[(all_sample_summary_cells$sample.ID %in% G12V_sections) & (all_sample_summary_cells$Variant.probe.name %in% G12V_key)])

data_1 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("G12D_cancer_cell", "non-G12D_cancer_and_Normal_cell"))

data_1["non-G12D_cancer_and_Normal_cell","WT"] <- normal_cell_REF_G12V_count
data_1["non-G12D_cancer_and_Normal_cell","ALT"] <- normal_cell_VAR_G12V_count
data_1["G12D_cancer_cell","WT"] <- G12D_cell_REF_G12V_count
data_1["G12D_cancer_cell","ALT"] <- G12D_cell_VAR_G12V_count
x <- data_1[,'ALT']
#print(x)
n <- rowSums(data_1)
#print(n)
variant_prop_test_1 <- prop.test(x, n, alternative = "greater", correct = FALSE) # 0.9452 G12D vs normal

data_2 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("G12V_cancer_cell", "non-G12D_cancer_and_Normal_cell"))
data_2["non-G12D_cancer_and_Normal_cell","WT"] <- normal_cell_REF_G12V_count
data_2["non-G12D_cancer_and_Normal_cell","ALT"] <- normal_cell_VAR_G12V_count
data_2["G12V_cancer_cell","WT"] <- G12V_cell_REF_G12V_count
data_2["G12V_cancer_cell","ALT"] <- G12V_cell_VAR_G12V_count
x <- data_2[,'ALT']
#print(x)
n <- rowSums(data_2)
#print(n)
variant_prop_test_2 <- prop.test(x, n, alternative = "greater", correct = FALSE) #2.2e-16 G12V vs normal

data_3 <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                      WT = c(0,0), # Cancer_cell, Normal_Cell
                      row.names = c("G12V_cancer_cell", "G12D_cancer_cell"))
data_3["G12D_cancer_cell","WT"] <- G12D_cell_REF_G12V_count
data_3["G12D_cancer_cell","ALT"] <- G12D_cell_VAR_G12V_count
data_3["G12V_cancer_cell","WT"] <- G12V_cell_REF_G12V_count
data_3["G12V_cancer_cell","ALT"] <- G12V_cell_VAR_G12V_count
x <- data_3[,'ALT']
#print(x)
n <- rowSums(data_3)
#print(n)
variant_prop_test_3 <- prop.test(x, n, alternative = "greater", correct = FALSE) # 2.2e-16 G12V vs G12D

normal_var_count <- normal_cell_VAR_G12V_count
normal_ref_count <- normal_cell_REF_G12V_count
clone_2_var_count <- G12V_cell_VAR_G12V_count
clone_2_ref_count <- G12V_cell_REF_G12V_count
clone_1_var_count <- G12D_cell_VAR_G12V_count
clone_1_ref_count <- G12D_cell_REF_G12V_count

normal_total = normal_var_count + normal_ref_count
clone_2_total = clone_2_var_count + clone_2_ref_count
clone_1_total = clone_1_var_count + clone_1_ref_count

n_cell_normal = normal_cell_VAR_G12V + normal_cell_REF_G12V
n_cell_clone_2 = G12V_cell_VAR_G12V + G12V_cell_REF_G12V
n_cell_clone_1 = G12D_cell_VAR_G12V + G12D_cell_REF_G12V

data_df <- data.frame(cell_type = rep(c("Normal cell", "G12V cancer cell", "G12D cancer cell"), each = 2),
                      proportion = c(normal_var_count/normal_total,normal_ref_count/normal_total,clone_2_var_count/clone_2_total,clone_2_ref_count/clone_2_total,clone_1_var_count/clone_1_total,clone_1_ref_count/clone_1_total),
                      prop_total = rep(c(1, 1, 1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "G12V cancer cell", "G12D cancer cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants) # 0.9452 G12D vs normal
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants) # 2.2e-16 G12V vs normal
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants) # 2.2e-16 G12V vs G12D
p1 <- ggplot(data_df, aes(x = cell_type, y = proportion, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1, by = 0.2)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "G12D cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.2, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("G12D cancer cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.4, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Proportion",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of KRAS p.G12V probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/KRAS_p.G12V_barplots_proportion_p_value_alt_g12d_vs_g12v.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

data_df <- data.frame(cell_type = rep(c("Normal cell", "G12V cancer cell", "G12D cancer cell"), each = 2),
                      probe_count = c(normal_var_count,normal_ref_count,clone_2_var_count,clone_2_ref_count,clone_1_var_count,clone_1_ref_count),
                      prop_total = rep(c(n_cell_normal,n_cell_clone_2,n_cell_clone_1), each = 2),
                      cell_count = rep(paste0("n_cell = ", c(n_cell_normal,n_cell_clone_2,n_cell_clone_1)), each = 2),
                      Allele_detected = rep(c("Variant", "Reference"), times = 3))
data_df$cell_type <- factor(data_df$cell_type, levels = c("Normal cell", "G12V cancer cell", "G12D cancer cell"))
data_df$Allele_detected <- factor(data_df$Allele_detected, levels = c("Variant", "Reference"))
number_of_variants = 42
prop_test_adj_p_c1_n <- p.adjust(variant_prop_test_1$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_n <- p.adjust(variant_prop_test_2$p.value, method = "bonferroni", n = number_of_variants)
prop_test_adj_p_c2_c1 <- p.adjust(variant_prop_test_3$p.value, method = "bonferroni", n = number_of_variants)
df_max_round = max(c(n_cell_normal,n_cell_clone_2,n_cell_clone_1))
p1 <- ggplot(data_df, aes(x = cell_type, y = probe_count, fill = Allele_detected)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(limits = c(0, 1.5*df_max_round)) +
    # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
    # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
    geom_signif(comparisons=list(c("Normal cell", "G12D cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c1_n, format = "E", digits = 2)),
                y_position = 1.2*df_max_round, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("G12D cancer cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_c1, format = "E", digits = 2)),
                y_position = 1.3*df_max_round, tip_length = 0, vjust=0) +
    geom_signif(comparisons=list(c("Normal cell", "G12V cancer cell")), 
                annotations=paste0("Adj prop test p = ", formatC(prop_test_adj_p_c2_n, format = "E", digits = 2)),
                y_position = 1.4*df_max_round, tip_length = 0, vjust=0) +
    geom_text(data = subset(data_df, Allele_detected == "Variant"),
              aes(label = cell_count, x = cell_type, y = prop_total),
              color = "black", nudge_y = max(data_df$prop_total)*0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 12),
          #legend.position="none",
          axis.title=element_text(size=12),
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
    labs( y = paste("Probe count",sep = ""),
          #x = "Treatment", 
          title = paste0("Proportion of KRAS p.G12V probes in cells")) + #,
    #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
    theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
pdf(paste0(out_dir,"/KRAS_p.G12V_barplots_count_p_value_alt_g12d_vs_g12v.pdf"),useDingbats = F, height=5, width = 5)
print(p1)
dev.off()

# Figure 6 DotPlot, hallmark enrichR pathways and associated module score VlnPlots
# Calculating differentially expressed genes in each cell type based on mutation status for the merged object.
# going to start by just looking at the results for tumor cells but I also want the results for other cell 
# types in case it makes a difference there.
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/"
prefix = "PDAC_merge_primary_KRAS_20241205_single_SCT"
meta.data <- read.table(paste0(out_dir,"/PDAC_merge_primary_KRAS_single_SCT_meta.data_20240102.tsv"),sep="\t",header=T)
obj <- readRDS(paste0(out_dir,"/",prefix,"_20240102.rds"))
dir.create(paste0(out_dir,"/","FindMarkers"))
print("cell_type_merge_v3")
cell_types <- unique(obj$cell_type_merge_v3)
print(cell_types)
obj$cell_type_merge_v3_KRAS_call <- paste0(obj$cell_type_merge_v3, "__", obj$Section_KRAS_call)
Idents(obj) <-"cell_type_merge_v3_KRAS_call"
for (cell_type in cell_types) {
    cell_type_genotypes <- unique(obj$cell_type_merge_v3_KRAS_call)[grepl(cell_type,unique(obj$cell_type_merge_v3_KRAS_call),fixed=T)]
    print(cell_type_genotypes)
    if (length(cell_type_genotypes) == 2) {
        print(cell_type)
        cell_type_G12D_label <- cell_type_genotypes[grepl("G12D",cell_type_genotypes)]
        cell_type_G12V_label <- cell_type_genotypes[grepl("G12V",cell_type_genotypes)]
        cell_type_degs <- FindMarkers(obj,assay = "SCT", slot="data", ident.1= cell_type_G12D_label, ident.2=cell_type_G12V_label, logfc.threshold = 0, min.pct=0.05, min.diff.pct=0.1, test.use='wilcox')
        write.table(cell_type_degs, paste0(out_dir,"/FindMarkers/",prefix,"_cell_type_merge_v3_",cell_type,"_id1_G12D_id2_G12V_20240106.tsv"),sep='\t',quote=F)
    }
}
# tumor has 64 degs
# Naive_T_cell has 24 degs
# every other cell type is either a very small number of degs or they are almost 
# all up or all down which is normally just a result of one specific sample 
print("cell_type_merge_v2")
cell_types <- unique(obj$cell_type_merge_v2)
print(cell_types)
obj$cell_type_merge_v2_KRAS_call <- paste0(obj$cell_type_merge_v2, "__", obj$Section_KRAS_call)
Idents(obj) <-"cell_type_merge_v2_KRAS_call"
for (cell_type in cell_types) {
    cell_type_genotypes <- unique(obj$cell_type_merge_v2_KRAS_call)[grepl(cell_type,unique(obj$cell_type_merge_v2_KRAS_call),fixed=T)]
    if (length(cell_type_genotypes) == 2) {
        print(cell_type)
        cell_type_G12D_label <- cell_type_genotypes[grepl("G12D",cell_type_genotypes)]
        cell_type_G12V_label <- cell_type_genotypes[grepl("G12V",cell_type_genotypes)]
        cell_type_degs <- FindMarkers(obj, assay = "SCT", slot="data", ident.1= cell_type_G12D_label, ident.2=cell_type_G12V_label, logfc.threshold = 0, min.pct=0.05, min.diff.pct=0.1, test.use='wilcox')
        write.table(cell_type_degs, paste0(out_dir, "/FindMarkers/", prefix, "_cell_type_merge_v2_", cell_type, "_id1_G12D_id2_G12V_20240106.tsv"), sep='\t', quote=F)
    }
}
### plotting difference in total cell type proportions across the first 100 layers moving away from tumor cells.
# Moffit subtyping
basal_geneset <- c('VGLL1','UCA1','S100A2','LY6D','SPRR3','SPRR1B','LEMD1','KRT15','CTSL2','DHRS9','AREG','CST6','SERPINB3','KRT6A','KRT6C',"SERPINB4",'FAM83A','SCEL','FGFBP1','KRT7','KRT17','GPR87','TNS4','SLC2A1','ANXA8L2')
classical_geneset <- c('BTNL8','FAM3D','ATAD4','AGR3','CTSE','LOC400573','LYZ','TFF2','TFF1','ANXA10','LGALS4','PLASG10','CEACAM6','VSIG2','TSPAN8','ST6GALNAC1','AGR2','TFF3','CYP3A7','MYO1A','CLRN3','KRT20','CDH17','SPINK4','REG4')
# conda activate banksy
sct_features <- Features(obj, assay = "SCT")
basal_geneset_xen <- c()
for (gene in basal_geneset) {
    if (gene %in% sct_features) {
        basal_geneset_xen <- c(basal_geneset_xen, gene)
    }
}
classical_geneset_xen <- c()
for (gene in classical_geneset) {
    if (gene %in% sct_features) {
        classical_geneset_xen <- c(classical_geneset_xen, gene)
    }
}
print(basal_geneset_xen)
# [1] "LY6D"     "SERPINB3" "FGFBP1"   "KRT7"    
print(classical_geneset_xen)
# [1] "AGR3"  "TFF2"  "KRT20"
obj <- AddModuleScore(obj, features = list(classical_geneset_xen), assay="SCT", name="classical_score", ctrl = 15, seed=1234) # you need to give AddModuleScore a list (e.g. list()) instead of a vector (e.g. c()) in order to get a single score for the entire set of genes. Otherwise, the function will generate one score for each gene.
classical_score_column <- colnames(obj@meta.data)[grepl("classical_score",colnames(obj@meta.data))]
obj <- AddModuleScore(obj, features = list(classical_geneset_xen), assay="SCT", name="basal_score", ctrl = 15, seed=1234) # you need to give AddModuleScore a list (e.g. list()) instead of a vector (e.g. c()) in order to get a single score for the entire set of genes. Otherwise, the function will generate one score for each gene.
basal_score_column <- colnames(obj@meta.data)[grepl("basal_score",colnames(obj@meta.data))]
module_score_columns <- c(classical_score_column,basal_score_column)
print(module_score_columns)
# "classical_score1" "basal_score1"
coordinates <- Embeddings(obj, reduction = "umap.50PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
DefaultAssay(obj) <- "SCT"
Tumor_g12_labels <- sort(unique(obj$cell_type_merge_v3_KRAS_call)[grepl("Tumor",unique(obj$cell_type_merge_v3_KRAS_call), fixed = TRUE)])
unified_cells_to_highlight <- list()
for (unified_subcluster_label in Tumor_g12_labels) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,"cell_type_merge_v3_KRAS_call"] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
}
print(length(unified_cells_to_highlight[[1]]))
# [1] 78461
print(length(unified_cells_to_highlight[[2]]))
# [1] 24176
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
DefaultAssay(obj) <- "SCT"
pdf(paste0(out_dir,"/",prefix,"_20240107_PDAC_subtype_module_score_dimplots.pdf"),useDingbats = F,height=5, width = 5)
print(rasterize(DimPlot(obj, group.by = "cell_type_merge_v3", reduction = "umap.50PC", pt.size=0.3, label=F, raster=FALSE)+ scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + labs(title = paste0("cell_type_merge_v3")), layers="Point", dpi=300))
print(rasterize(DimPlot(obj, group.by = "cell_type_merge_v3", reduction = "umap.50PC", pt.size=0.3, label=F, raster=FALSE)+ scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300))
for (i in 1:length(Tumor_g12_labels)) {
    unified_subcluster_label <- Tumor_g12_labels[i]
    print(unified_subcluster_label)
    cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = "cell_type_merge_v3_KRAS_call", cells.highlight = cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=F, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = "cell_type_merge_v3_KRAS_call", cells.highlight = cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=F, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers='Point', dpi=300))
}
for (i in 1:length(module_score_columns)) {
    print(rasterize(FeaturePlot(obj, slot="counts", features = module_score_columns[i], min.cutoff= 0, max.cutoff="q95", order = T, reduction = "umap.50PC", label=F, pt.size=0.3, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300))
    print(rasterize(FeaturePlot(obj, slot="counts", features = module_score_columns[i], min.cutoff= 0, max.cutoff="q95", order = T, reduction = "umap.50PC", label=F, pt.size=0.3, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300))
}
KRAS_g12_counts <- c("KRAS.p.G12D.ALT.T_Xenium.with.snvs_count","KRAS.p.G12V.ALT.A_Xenium.with.snvs_count","KRAS.p.G12V.WT_Xenium.with.snvs_count")
for (i in 1:length(module_score_columns)) {
    print(rasterize(FeaturePlot(obj, slot="counts", features = KRAS_g12_counts[i], min.cutoff= 0, max.cutoff="q95", order = T, reduction = "umap.50PC", label=F, pt.size=0.3, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300))
    print(rasterize(FeaturePlot(obj, slot="counts", features = KRAS_g12_counts[i], min.cutoff= 0, max.cutoff="q95", order = T, reduction = "umap.50PC", label=F, pt.size=0.3, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
}
dev.off()
DefaultAssay(obj) <- "SCT"
Idents(obj) <- "cell_type_merge_v3"
pdf(paste0(out_dir,"/",prefix,"_20240107_PDAC_subtype_module_score_vlnplots.pdf"),useDingbats = F,height=5, width = 5)
for (i in 1:length(module_score_columns)) {
    print(VlnPlot(obj, features = module_score_columns[i], assay="SCT", split.by = "Section_KRAS_call", idents = c("Tumor"), pt.size=0, raster=FALSE)) #, max.cutoff="q10")# first generate un-rastered plot or the "order = T" argument is ignored
    #print(VlnPlot(obj, features = module_score_columns[i], assay="SCT", split.by = "Section_KRAS_call", idents = c("Tumor"), pt.size=0, raster=FALSE))
}
dev.off()
hallmark_tables <- read.table("/diskmnt/Projects/Users/austins2/tools/gene_tables/MSigDB/HALLMARK/all_HALLMARK_geneset_tables.tsv", sep='\t',header=T, row.names = "index")
hallmark_gene_lists <- list()
geneset_names <- hallmark_tables$geneset_name
for (pathway_name in geneset_names) {
    geneset <- readLines(hallmark_tables[(hallmark_tables$geneset_name == pathway_name),"geneset_path"])
    hallmark_gene_lists[[pathway_name]] <- geneset
}
meta.data.checkpoint <- obj@meta.data
sct_features <- Features(obj, assay = "SCT")
Tumor_degs <- read.table(paste0(out_dir,"/FindMarkers/",prefix,"_cell_type_merge_v3_",'Tumor',"_id1_G12D_id2_G12V_20240106.tsv"),sep='\t',header=T)
Tumor_degs$gene <- rownames(Tumor_degs)
Tumor_degs$pathway <- NA
for (gene in rownames(Tumor_degs)) {
    for (pathway_name in geneset_names) {
        geneset <- hallmark_gene_lists[[pathway_name]]
        if (gene %in% geneset) {
            if (is.na(Tumor_degs[gene,"pathway"])) {
                Tumor_degs[gene,"pathway"] <- pathway_name
            } else {
                Tumor_degs[gene,"pathway"] <- paste(c(noquote(pathway_name),noquote(Tumor_degs[gene,"pathway"])), collapse = ',')
            }
        }
    }
}
write.table(Tumor_degs,paste0(out_dir,"/FindMarkers/",prefix,"_cell_type_merge_v3_",'Tumor',"_id1_G12D_id2_G12V_20240106_with_hallmark_pathways.tsv"),sep='\t',quote=F)
cell_types_kras <- unique(obj$cell_type_merge_v3_KRAS_call)
mean_module_scores_mat <- matrix(0, length(geneset_names), length(cell_types_kras))
rownames(mean_module_scores_mat) <- geneset_names
colnames(mean_module_scores_mat) <- cell_types_kras
geneset_to_include_for_plotting <- c()
geneset_names_with_gene_count <- c()
for (pathway_name in geneset_names) {
    print(pathway_name)
    geneset <- hallmark_gene_lists[[pathway_name]]
    geneset_xen <- c()
    for (gene in geneset) {
        if (gene %in% sct_features) {
            geneset_xen <- c(geneset_xen, gene)
        }
    }
    print(length(geneset_xen))
    print(geneset_xen)
    geneset_names_with_gene_count <- c(geneset_names_with_gene_count, paste0(pathway_name,"__",length(geneset_xen)))
    if (length(geneset_xen) > 0) {
        if (length(geneset_xen) > 2) {
            geneset_to_include_for_plotting <- c(geneset_to_include_for_plotting, pathway_name)
        }
        obj <- AddModuleScore(obj, features = list(geneset_xen), assay="SCT", name=pathway_name, ctrl = 15, seed=1234)
        for (cell_type_kras in cell_types_kras) {
            metadata_columns <- colnames(obj@meta.data)
            module_name_meta <- metadata_columns[grepl(pathway_name, metadata_columns)]
            cell_kras_module_scores = obj@meta.data[(obj$cell_type_merge_v3_KRAS_call == cell_type_kras),module_name_meta]
            cell_kras_module_mean = mean(cell_kras_module_scores)
            mean_module_scores_mat[rownames(mean_module_scores_mat) == pathway_name, colnames(mean_module_scores_mat) == cell_type_kras] <- cell_kras_module_mean
        }
    }
}
# [1] "HALLMARK_ADIPOGENESIS"                                                                                                                                              
# [1] 9                                                                                                                                                                    
# [1] "ADIPOQ" "CYP4B1" "LPL"    "MYLK"   "PPARG"  "CAVIN1" "RETN"   "CAVIN2"                                                                                              
# [9] "SNCG"                                                                                                                                                               
# [1] "HALLMARK_ALLOGRAFT_REJECTION"                                                                                                                                       
# [1] 27                                                                                                                                                                   
# [1] "CCL19" "CCL5"  "CCR2"  "CD2"   "CD247" "CD28"  "CD3D"  "CD3E"  "CD4"                                                                                               
# [10] "CD79A" "CD86"  "CD8A"  "CXCL9" "EGFR"  "FAS"   "GZMA"  "GZMB"  "IGSF6"                                                                                             
# [19] "IL2RA" "IRF8"  "KLRD1" "LIF"   "LY86"  "PRF1"  "PTPRC" "SPI1"  "THY1"
# [1] "HALLMARK_ANDROGEN_RESPONSE"                                                                                                                                         
# [1] 5                                                                                                                                                                    
# [1] "ADAMTS1" "ALDH1A3" "SLC26A2" "SPDEF"   "STEAP4"                                                                                                                     
# [1] "HALLMARK_ANGIOGENESIS"                                                                                                                                              
# [1] 5                                                                                                                                                                    
# [1] "COL5A2" "CXCL6"  "LPL"    "STC1"   "VCAN"                                                                                                                           
# [1] "HALLMARK_APICAL_JUNCTION"                                                                                                                                           
# [1] 12                                                                                                                                                                   
# [1] "ACTG2"   "CD274"   "CD34"    "CD86"    "COL17A1" "EGFR"    "FBN1"                                                                                                  
# [8] "PECAM1"  "PTPRC"   "THY1"    "VCAN"    "VWF"                                                                                                                       
# [1] "HALLMARK_APICAL_SURFACE"                                                                                                                                            
# [1] 3                                                                                                                                                                    
# [1] "GHRL" "SRPX" "THY1"                                                                                                                                                 
# [1] "HALLMARK_APOPTOSIS"                                                                                                                                                 
# [1] 10                                                                                                                                                                   
# [1] "BCL2L11" "CAV1"    "CD14"    "CD2"     "CD69"    "ERBB2"   "FAS"                                                                                                   
# [8] "PDGFRB"  "PRF1"    "TOP2A"                                                                                                                                         
# [1] "HALLMARK_BILE_ACID_METABOLISM"                                                                                                                                      
# [1] 3                                                                                                                                                                    
# [1] "AQP9"  "AR"    "BBOX1" 
# [1] "HALLMARK_CHOLESTEROL_HOMEOSTASIS"                                                                                                                                   
# [1] 3                                                                                                                                                                    
# [1] "ADH4"  "LPL"   "PPARG"                                                                                                                                              
# [1] "HALLMARK_COAGULATION"                                                                                                                                               
# [1] 9                                                                                                                                                                    
# [1] "CFB"      "CTSK"     "FBN1"     "HMGCS2"   "PECAM1"   "RAPGEF3"  "S100A1"                                                                                           
# [8] "SERPINB2" "VWF"                                                                                                                                                     
# [1] "HALLMARK_COMPLEMENT"                                                                                                                                                
# [1] 9                                                                                                                                                                    
# [1] "CCL5"     "CFB"      "FCN1"     "GZMA"     "GZMB"     "GZMK"     "PLA2G7"                                                                                           
# [8] "S100A12"  "SERPINB2"                                                                                                                                                
# [1] "HALLMARK_DNA_REPAIR"                                                                                                                                                
# [1] 1                                                                                                                                                                    
# [1] "PCNA"                                                                                                                                                               
# [1] "HALLMARK_E2F_TARGETS"                                                                                                                                               
# [1] 6                                                                                                                                                                    
# [1] "CCNB2" "CDK1"  "MKI67" "MYC"   "PCNA"  "TOP2A"                                                                                                                      
# [1] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"                                                                                                                         
# [1] 24                                                                                                                                                                   
# [1] "ACTA2"  "ANPEP"  "BASP1"  "COL5A2" "CXCL6"  "DST"    "FAS"    "FBLN1"                                                                                              
# [9] "FBN1"   "FSTL3"  "GEM"    "GLIPR1" "GPC1"   "MEST"   "MFAP5"  "MYLK"                                                                                               
# [17] "PCOLCE" "PDGFRB" "PMP22"  "SFRP4"  "THBS2"  "THY1"   "TNC"    "VCAN"
# [1] "HALLMARK_ESTROGEN_RESPONSE_EARLY"                                                                                                                                   
# [1] 8                                                                                                                                                                    
# [1] "AQP3"    "AR"      "FHL2"    "MLPH"    "MYC"     "PGR"     "SLC26A2"                                                                                                
# [8] "STC2"                                                                                                                                                               
# [1] "HALLMARK_ESTROGEN_RESPONSE_LATE"                                                                                                                                    
# [1] 10                                                                                                                                                                   
# [1] "ASCL1"   "CAV1"    "DUSP2"   "HMGCS2"  "KLK11"   "MEST"    "PCP4"                                                                                                  
# [8] "PGR"     "SLC26A2" "TOP2A"                                                                                                                                         
# [1] "HALLMARK_FATTY_ACID_METABOLISM"                                                                                                                                     
# [1] 5                                                                                                                                                                    
# [1] "ADH1C"  "CA4"    "CYP1A1" "HMGCS2" "INMT"                                                                                                                           
# [1] "HALLMARK_G2M_CHECKPOINT"                                                                                                                                            
# [1] 7                                                                                                                                                                    
# [1] "CCNB2" "CDK1"  "CENPF" "MKI67" "MYC"   "TOP2A" "UBE2C"                                                                                                              
# [1] "HALLMARK_GLYCOLYSIS"                                                                                                                                                
# [1] 9                                                                                                                                                                    
# [1] "CDK1"  "CXCR4" "EGFR"  "GPC1"  "GPC3"  "MET"   "STC1"  "STC2"  "VCAN"                                                                                               
# [1] "HALLMARK_HEDGEHOG_SIGNALING"                                                                                                                                        
# [1] 1                                                                                                                                                                    
# [1] "THY1"
# [1] "HALLMARK_HEME_METABOLISM"                                                                                                                                           
# [1] 9                                                                                                                                                                    
# [1] "AHSP"   "ALAS2"  "AQP3"   "ACKR1"  "TENT5C" "GYPA"   "GYPB"   "SLC4A1"                                                                                              
# [9] "SNCA"                                                                                                                                                               
# [1] "HALLMARK_HYPOXIA"                                                                                                                                                   
# [1] 9                                                                                                                                                                    
# [1] "CAV1"   "CXCR4"  "EGFR"   "GPC1"   "GPC3"   "CAVIN1" "SRPX"   "STC1"                                                                                                
# [9] "STC2"                                                                                                                                                               
# [1] "HALLMARK_IL2_STAT5_SIGNALING"                                                                                                                                       
# [1] 20                                                                                                                                                                   
# [1] "AGER"    "CD83"    "CD86"    "COCH"    "CTLA4"   "CXCL10"  "EOMES"                                                                                                 
# [8] "FGL2"    "IL1R2"   "IL1RL1"  "IL2RA"   "IL3RA"   "IRF8"    "ITGAE"                                                                                                 
# [15] "LIF"     "MYC"     "RGS16"   "S100A1"  "SELL"    "TNFRSF9"                                                                                                         
# [1] "HALLMARK_IL6_JAK_STAT3_SIGNALING"                                                                                                                                   
# [1] 8                                                                                                                                                                    
# [1] "CD14"   "CSF2RA" "CXCL10" "CXCL9"  "FAS"    "IL1R2"  "IL2RA"  "IL3RA"                                                                                               
# [1] "HALLMARK_INFLAMMATORY_RESPONSE"                                                                                                                                     
# [1] 26                                                                                                                                                                   
# [1] "AQP9"    "CCL5"    "CCR7"    "CD14"    "CD69"    "CD70"    "CSF3"                                                                                                  
# [8] "CXCL10"  "CXCL6"   "CXCL9"   "EDN1"    "ADGRE1"  "GPC3"    "GPR183"                                                                                                
# [15] "IL7R"    "LAMP3"   "LIF"     "MARCO"   "MET"     "MYC"     "PDPN"                                                                                                  
# [22] "RGS16"   "SELE"    "SELL"    "SLAMF1"  "TNFRSF9"
# [1] "HALLMARK_INTERFERON_ALPHA_RESPONSE"                                                                                                                                 
# [1] 3                                                                                                                                                                    
# [1] "CXCL10" "LAMP3"  "SELL"                                                                                                                                             
# [1] "HALLMARK_INTERFERON_GAMMA_RESPONSE"                                                                                                                                 
# [1] 14                                                                                                                                                                   
# [1] "BANK1"  "CCL5"   "CD274"  "CD69"   "CD86"   "CFB"    "CXCL10" "CXCL9"                                                                                              
# [9] "FAS"    "FCGR1A" "FGL2"   "GZMA"   "IRF8"   "SLAMF7"                                                                                                               
# [1] "HALLMARK_KRAS_SIGNALING_DN"                                                                                                                                         
# [1] 8                                                                                                                                                                    
# [1] "CDH16"    "EDN1"     "TENT5C"   "PDCD1"    "SERPINB2" "TCL1A"    "TFF2"                                                                                             
# [8] "UPK3B"                                                                                                                                                              
# [1] "HALLMARK_KRAS_SIGNALING_UP"                                                                                                                                         
# [1] 19                                                                                                                                                                   
# [1] "ALDH1A3" "CFB"     "CSF2RA"  "CXCL10"  "CXCR4"   "ADGRL4"  "GNG11"                                                                                                 
# [8] "IL7R"    "IRF8"    "LIF"     "MALL"    "PCP4"    "PECAM1"  "PRDM1"                                                                                                 
# [15] "RETN"    "RGS16"   "TFPI"    "TMEM100" "VWA5A"                                                                                                                     
# [1] "HALLMARK_MITOTIC_SPINDLE"                                                                                                                                           
# [1] 6                                                                                                                                                                    
# [1] "BCL2L11" "CCNB2"   "CDK1"    "CENPF"   "DST"     "TOP2A"
# [1] "HALLMARK_MTORC1_SIGNALING"                                                                                                                                          
# [1] 3                                                                                                                                                                    
# [1] "CXCR4" "FGL2"  "STC1"                                                                                                                                               
# [1] "HALLMARK_MYC_TARGETS_V1"                                                                                                                                            
# [1] 2                                                                                                                                                                    
# [1] "MYC"  "PCNA"                                                                                                                                                        
# [1] "HALLMARK_MYC_TARGETS_V2"                                                                                                                                            
# [1] 2                                                                                                                                                                    
# [1] "DUSP2" "MYC"                                                                                                                                                        
# [1] "HALLMARK_MYOGENESIS"                                                                                                                                                
# [1] 8                                                                                                                                                                    
# [1] "DES"   "IGF1"  "MEF2C" "MYH11" "MYLK"  "PVALB" "SPDEF" "STC2"                                                                                                       
# [1] "HALLMARK_NOTCH_SIGNALING"                                                                                                                                           
# [1] 0                                                                                                                                                                    
# NULL                                                                                                                                                                     
# [1] "HALLMARK_OXIDATIVE_PHOSPHORYLATION"                                                                                                                                 
# [1] 0                                                                                                                                                                    
# NULL                                                                                                                                                                     
# [1] "HALLMARK_P53_PATHWAY"                                                                                                                                               
# [1] 8                                                                                                                                                                    
# [1] "CLCA2" "FAS"   "GPX2"  "LIF"   "MDM2"  "PCNA"  "RGS16" "VWA5A"
# [1] "HALLMARK_PANCREAS_BETA_CELLS"                                                                                                                                       
# [1] 6                                                                                                                                                                    
# [1] "CHGA"  "GCG"   "INS"   "PCSK2" "SCGN"  "SST"                                                                                                                        
# [1] "HALLMARK_PEROXISOME"                                                                                                                                                
# [1] 2                                                                                                                                                                    
# [1] "SEMA3C" "TOP2A"                                                                                                                                                     
# [1] "HALLMARK_PI3K_AKT_MTOR_SIGNALING"                                                                                                                                   
# [1] 3                                                                                                                                                                    
# [1] "CDK1"  "CXCR4" "EGFR"                                                                                                                                               
# [1] "HALLMARK_PROTEIN_SECRETION"                                                                                                                                         
# [1] 2                                                                                                                                                                    
# [1] "DST"  "EGFR"                                                                                                                                                        
# [1] "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"                                                                                                                           
# [1] 1                                                                                                                                                                    
# [1] "MBP"                                                                                                                                                                
# [1] "HALLMARK_SPERMATOGENESIS"                                                                                                                                           
# [1] 3                                                                                                                                                                    
# [1] "CCNB2" "CDK1"  "CFTR"                                                                                                                                               
# [1] "HALLMARK_TGF_BETA_SIGNALING"                                                                                                                                        
# [1] 1                                                                                                                                                                    
# [1] "LTBP2"
# [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# [1] 16
# [1] "CCL5"     "CD69"     "CD83"     "CXCL10"   "CXCL2"    "CXCL6"   
# [7] "DUSP2"    "EDN1"     "GEM"      "GPR183"   "IL7R"     "LIF"     
# [13] "MYC"      "SERPINB2" "TNC"      "TNFRSF9" 
# [1] "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
# [1] 1
# [1] "STC2"
# [1] "HALLMARK_UV_RESPONSE_DN"
# [1] 12
# [1] "CAV1"   "COL5A2" "ERBB2"  "FHL2"   "KCNMA1" "KIT"    "MET"    "MYC"   
# [9] "PDGFRB" "PMP22"  "PPARG"  "TFPI"  
# [1] "HALLMARK_UV_RESPONSE_UP"
# [1] 5
# [1] "AQP3"    "BCL2L11" "CXCL2"   "CYP1A1"  "EPCAM"  
# [1] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
# [1] 1
# [1] "MYC"
# [1] "HALLMARK_XENOBIOTIC_METABOLISM"
# [1] 11
# [1] "ADH1C"  "AQP9"   "CFB"    "CYP1A1" "ESR1"   "FAS"    "FBLN1"  "IGF1"  
# [9] "IRF8"   "PTGDS"  "TAT" 
rownames(mean_module_scores_mat) <- geneset_names_with_gene_count
dir.create(paste0(out_dir,"/","KRAS_var_pathway_analysis","/"))
write.table(mean_module_scores_mat, paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_mean_pathway_module_scores.tsv"),sep='\t',quote=F)
# based on the results in this table:
# # HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION - higher in G12D (pct.1)
# TNC
# PMP22
# ANPEP
# COL5A2
# VCAN
# # HALLMARK_UV_RESPONSE_DN - higher in G12D (pct.1)
# CAV1
# TFPI
# MYC
# MET
# PMP22
# COL5A2
# # # HALLMARK_HYPOXIA - higher in G12D (pct.1)
# # CAVIN1
# # CAV1
# # SRPX
# # HALLMARK_TNFA_SIGNALING_VIA_NFKB - higher in G12V (pct.2)
# LIF
# CXCL2
# CXCL6
# EDN1
# SERPINB2
# # HALLMARK_P53_PATHWAY
# LIF - negative regulator of TP53: PMID: 26161442
# GPX2 - considered to be in the P53 pathway but is actually regulated by p63 (expression is induced by p63 activity), however "Activation of GPX2 gene can alleviate hydrogen peroxide-induced cell death via p53 activation and it also indicates that GPX2 activation can promote cancer cell growth."
# PCNA - "suggest a complex cellular response to DNA damage in which p53 transiently activates expression of PCNA for the purpose of limited DNA repair." https://doi.org/10.1128/mcb.19.1.12 and "demonstrate p53 binding to a target site in the PCNA promoter, recruitment of p300/CREB-binding protein, and localized acetylation of histone H4 in an IR-dependent manner. These molecular events are likely to play a role in mediating activation of PCNA gene expression by p53 during the cellular response to DNA damage." - https://doi.org/10.1074/jbc.M302671200
# FAS -  "Earlier work established that p53 target promoters display pronounced differences in RNAPII occupancy prior to p53 activation, with higher levels observed at cell cycle arrest genes (e.g., CDKN1A) relative to apoptotic genes (e.g., FAS and BBC3), which correlated with a delayed in the accumulation of mature FAS mRNAs."- https://doi.org/10.1038/cdd.2017.174 and "p53 induces the expression of the death receptors Fas/Fas ligand and KILLER/DR5 located on the cell membrane, which activates caspase 8 and leads to apoptosis" - https://doi.org/10.1038/s41392-023-01347-1
# MDM2 "E3 ubiquitin-protein ligase that mediates ubiquitination of p53/TP53, leading to its degradation by the proteasome" <- UniProtKB
# VWA5A
# SMYD2 <- Reactome Pathways 2024: "Regulation of TP53 Activity Through Methylation" # also described in UniProtKB as suppressing TP53 activity via histone methylation
# HALLMARK_G2M_CHECKPOINT
# "CCNB2" 
# "CDK1"
# "CENPF"
# "MKI67"
# "MYC"
# "TOP2A"
# "UBE2C"    
# # HALLMARK_INFLAMMATORY_RESPONSE
# SELL - lymphocyte adherence to endothelial - NCBI Gene Summary
# MYC 
# LAMP3 - macrophages and dendritic cells
# RGS16 
# GPR183 
# CCR7 
# CD14 
# IL7R 
# MET 
# make a violin plot of geneset expression in just the tumor cells
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_HALLMARK_pathway_module_scores.pdf"), useDingbats = F,height=5, width = 5)
Idents(obj) <- "cell_type_merge_v3"
for (pathway_name in geneset_to_include_for_plotting) {
    print(pathway_name)
    p1 <- VlnPlot(obj, features = paste0(pathway_name,"1"), split.by = "Section_KRAS_call", idents = c("Tumor"), pt.size=0, raster=FALSE)
    scores_g12d <- obj@meta.data[(obj$Section_KRAS_call == "KRAS-p-G12D-ALT-T"), paste0(pathway_name,"1")]
    scores_g12v <- obj@meta.data[(obj$Section_KRAS_call == "KRAS-p-G12V-ALT-A"), paste0(pathway_name,"1")]
    w.test <- wilcox.test(scores_g12d, scores_g12v, alternative = "two.sided", exact = TRUE)
    exact_pvalue <- format.pval(w.test$p.value, digits = 3)
    print(exact_pvalue)
    p1 <- p1 + labs(subtitle= paste0("Wilcoxon p = ",exact_pvalue)) +
        theme(plot.title=element_text(size=10, hjust=0.5, face="bold", colour="black", vjust=1)) +
        theme(plot.subtitle=element_text(size=10, hjust=0.5, face="italic", color="black"))
    print(p1)
}
dev.off()
combined_geneset_xen <- c()
for (pathway_name in geneset_names) {
    print(pathway_name)
    geneset <- hallmark_gene_lists[[pathway_name]]
    geneset_xen <- c()
    for (gene in geneset) {
        if (gene %in% sct_features) {
            geneset_xen <- c(geneset_xen, gene)
        }
    }
    combined_geneset_xen <- c(combined_geneset_xen, geneset_xen)
}
hallmark_genes <- AggregateExpression(obj, assays = "SCT", features = combined_geneset_xen, group.by=group_ident)
hallmark_genes_v2 = as.data.frame(aggregate_genes_result$SCT)
agg_hall_rows <- rownames(hallmark_genes_v2)
tumor_hallmark_rows <- agg_hall_rows[grepl("Tumor",agg_hall_rows)]
library(enrichR)
Idents(obj) <- "cell_type_merge_v3_KRAS_call"
test <- DEenrichRPlot(obj, ident.1 = "Tumor__KRAS-p-G12D-ALT-T", ident.2 = "Tumor__KRAS-p-G12V-ALT-A", 
                      assay = "SCT", test.use = "wilcox", enrich.database = "MSigDB_Hallmark_2020", balanced = T, 
                      max.genes = 100, cols = c("blue", "white", "red"), return.gene.list = TRUE) # this confirms the deg mapping to HALLMARK pathways above
write.table(test$pos, paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_HALLMARK_DEenrichRPlot_pos_enriched_in_G12D.tsv"),sep='\t',quote=F)
write.table(test$neg, paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_HALLMARK_DEenrichRPlot_neg_enriched_in_G12V.tsv"),sep='\t',quote=F)
p1 <- DEenrichRPlot(obj, ident.1 = "Tumor__KRAS-p-G12D-ALT-T", ident.2 = "Tumor__KRAS-p-G12V-ALT-A", # can also visualize stuff here: https://maayanlab.cloud/Enrichr/
                    assay = "SCT", test.use = "wilcox", enrich.database = "MSigDB_Hallmark_2020", balanced = T, 
                    max.genes = 100, cols = c("blue", "white", "red")) 
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_HALLMARK_DEenrichRPlot.pdf"), useDingbats = F,height=5, width = 15)
print(p1) # include this in the supplement. This 
dev.off()
deg.er.react.g12 <- DEenrichRPlot(obj, ident.1 = "Tumor__KRAS-p-G12D-ALT-T", ident.2 = "Tumor__KRAS-p-G12V-ALT-A", 
                                  assay = "SCT", test.use = "wilcox", enrich.database = "Reactome_Pathways_2024", balanced = T,  # can also visualize stuff here: https://maayanlab.cloud/Enrichr/
                                  max.genes = 100, cols = c("blue", "white", "red"), return.gene.list = TRUE) # this suggest that there likely is more negative regulation of TP53 in the p.G12V sections (Regulation of TP53 Activity Through Methylation)
write.table(deg.er.react.g12$pos, paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_Reactome_Pathways_2024_DEenrichRPlot_pos_enriched_in_G12D.tsv"),sep='\t',quote=F)
write.table(deg.er.react.g12$neg, paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_Reactome_Pathways_2024_DEenrichRPlot_neg_enriched_in_G12V.tsv"),sep='\t',quote=F)
p2 <- DEenrichRPlot(obj, ident.1 = "Tumor__KRAS-p-G12D-ALT-T", ident.2 = "Tumor__KRAS-p-G12V-ALT-A", 
                    assay = "SCT", test.use = "wilcox", enrich.database = "Reactome_Pathways_2024 ", balanced = T, 
                    max.genes = 100, cols = c("blue", "white", "red")) 
# didn't need to rely on the solution in this post but if it ends up happening that the neg and pos results are the same then you can follow here to ensure the positive and negative results are not identical.
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_Reactome_Pathways_2024_DEenrichRPlot.pdf"), useDingbats = F,height=5, width = 15)
print(p2) # include this in the supplement. This 
dev.off()
library(ggrastr)
library(ggpubr)
library(ggh4x)
genes_to_plot <- c('TNC','PMP22','ANPEP','COL5A2','VCAN','CAV1','TFPI','MYC','MET','PMP22','COL5A2','CAVIN1','CAV1',"SELL",'MYC','LAMP3','RGS16','GPR183','CCR7','CD14','IL7R','MET','LIF','CXCL2','CXCL6','EDN1','SERPINB2','LIF','GPX2','PCNA','FAS','MDM2','VWA5A',"CCNB2","CDK1","CENPF","MKI67","TOP2A","UBE2C")
genes_to_plot_uni <- unique(genes_to_plot)
print(genes_to_plot_uni)
#  [1] "TNC"      "PMP22"    "ANPEP"    "COL5A2"   "VCAN"     "CAV1"    
#  [7] "TFPI"     "MYC"      "MET"      "CAVIN1"   "LIF"      "CXCL2"   
# [13] "CXCL6"    "EDN1"     "SERPINB2" "GPX2"     "PCNA"     "FAS"     
# [19] "MDM2"     "VWA5A"    "CCNB2"    "CDK1"     "CENPF"    "MKI67"   
# [25] "TOP2A"    "UBE2C"    "SELL"     "LAMP3"    "RGS16"    "GPR183"  
# [31] "CCR7"     "CD14"     "IL7R"
obj_tumor_only <- subset(obj, subset = cell_type_merge_v3 == "Tumor")
obj_tumor_only <- SCTransform(obj_tumor_only, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2", conserve.memory = F)
p2 <- DotPlot(obj_tumor_only, assay = "SCT", features = genes_to_plot_uni, col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "Section_KRAS_call") + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red") + RotatedAxis() +
    force_panelsizes(rows = unit(0.4, "in"),
                     cols = unit(5.6, "in"))
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_DotPlot_MSigDB_Hallmark_2020.pdf"), useDingbats = F, height=4, width = 11)
print(p2)
dev.off()
Idents(obj) <- "cell_type_merge_v3_KRAS_call"
p3 <- DotPlot(obj, assay = "SCT", features = genes_to_plot_uni, col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "cell_type_merge_v3_KRAS_call", idents = c("Tumor__KRAS-p-G12V-ALT-A","Tumor__KRAS-p-G12D-ALT-T")) + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red") + RotatedAxis() +
    force_panelsizes(rows = unit(0.4, "in"),
                     cols = unit(5.6, "in"))
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_DotPlot_full_object_MSigDB_Hallmark_2020.pdf"), useDingbats = F, height=4, width = 11)
print(p3)
dev.off()
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_HALLMARK_pathway_module_scores_tumor_only_subset_renormalized.pdf"), useDingbats = F,height=5, width = 5)
Idents(obj) <- "cell_type_merge_v3"
for (pathway_name in geneset_to_include_for_plotting) {
    print(pathway_name)
    geneset <- hallmark_gene_lists[[pathway_name]]
    geneset_xen <- c()
    for (gene in geneset) {
        if (gene %in% sct_features) {
            geneset_xen <- c(geneset_xen, gene)
        }
    }
    obj <- AddModuleScore(obj, features = list(geneset_xen), assay="SCT", name=pathway_name, ctrl = 15, seed=1234)
    p1 <- VlnPlot(obj, features = paste0(pathway_name,"1"), split.by = "Section_KRAS_call", idents = c("Tumor"), pt.size=0, raster=FALSE)
    scores_g12d <- obj@meta.data[(obj$Section_KRAS_call == "KRAS-p-G12D-ALT-T"), paste0(pathway_name,"1")]
    scores_g12v <- obj@meta.data[(obj$Section_KRAS_call == "KRAS-p-G12V-ALT-A"), paste0(pathway_name,"1")]
    w.test <- wilcox.test(scores_g12d, scores_g12v, alternative = "two.sided", exact = TRUE)
    exact_pvalue <- format.pval(w.test$p.value, digits = 3)
    print(exact_pvalue)
    p1 <- p1 + labs(subtitle= paste0("Wilcoxon p = ",exact_pvalue)) +
        theme(plot.title=element_text(size=10, hjust=0.5, face="bold", colour="black", vjust=1)) +
        theme(plot.subtitle=element_text(size=10, hjust=0.5, face="italic", color="black"))
    print(p1)
}
dev.off()
write.table(obj@meta.data, paste0(out_dir,"/",prefix,"meta.data_20240113_pathways.tsv"),sep='\t',quote=F)
saveRDS(obj_tumor_only, paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_tumor_only_subset_20240113.rds"))
genes_to_plot_uni = c('TNC','PMP22','ANPEP','COL5A2','VCAN','CAV1','TFPI','MET','CAVIN1','LAMP3','CD14','CXCL6','EDN1','VWA5A')

p3 <- DotPlot(obj, assay = "SCT", features = genes_to_plot_uni, cols = "RdBu", col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "cell_type_merge_v3_KRAS_call", idents = c("Tumor__KRAS-p-G12V-ALT-A","Tumor__KRAS-p-G12D-ALT-T")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    force_panelsizes(rows = unit(0.4, "in"),
                     cols = unit(2.3, "in"))
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_DotPlot_full_object_MSigDB_Hallmark_2020_scRNA_validated.pdf"), useDingbats = F, height=4, width = 11)
print(p3)
dev.off()
p3 <- DotPlot(obj, assay = "SCT", features = genes_to_plot_uni, cols = "RdBu", col.min = -1, col.max = 1, dot.min = 0, dot.scale = 6, scale.min = 0, scale.max= 100, group.by = "cell_type_merge_v3_KRAS_call", split.by = "sample_ID", idents = c("Tumor__KRAS-p-G12V-ALT-A","Tumor__KRAS-p-G12D-ALT-T")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    force_panelsizes(rows = unit(3, "in"),
                     cols = unit(2.3, "in"))
pdf(paste0(out_dir,"/","KRAS_var_pathway_analysis","/",prefix,"_DotPlot_full_object_MSigDB_Hallmark_2020_scRNA_validated_sampleID.pdf"), useDingbats = F, height=4, width = 11)
print(p3)
dev.off()

#plotting the merged umap
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/"
obj_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/"
prefix = "PDAC_merge_primary_KRAS_20241205_single_SCT"
meta.data <- read.table(paste0(obj_dir,"/",prefix,"_meta.data_20240113_pathway_and_morph_tme_layers.tsv"),sep="\t",header=T)
obj <- readRDS("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only//PDAC_merge_primary_KRAS_20241205_single_SCT_20240102.rds")
obj$cell_type_merge_v3_reduced <- obj$cell_type_merge_v3
obj$cell_type_merge_v3_reduced[obj$cell_type_merge_v3 %in% c('CD4+T_cell','CD8+T_cell','Treg','Naive_T_cell')] <- "T_cell"
obj$cell_type_merge_v3_reduced[obj$cell_type_merge_v3 %in% c('myCAF','iCAF')] <- "Fibroblast"
obj$cell_type_merge_v3_reduced[obj$cell_type_merge_v3 %in% c('Islet_alpha','Islet_beta','Islet_delta','Islet_gamma','Islet_INS+GCG+')] <- "Islet"
obj$cell_type_merge_v3_reduced[obj$cell_type_merge_v3 %in% c('Ductal_reactive','Duct_like_1','Duct_like_2')] <- "Ductal"
# cell_type_group[["T_cell"]] <- c('CD4+T_cell','CD8+T_cell','Treg','Naive_T_cell')
# cell_type_group[["Fibroblast"]] <- c('myCAF','iCAF')
# cell_type_group[["Islet"]] <- c('Islet_alpha','Islet_beta','Islet_delta','Islet_gamma','Islet_INS+GCG+')
# cell_type_group[["Epithelial"]] <- c('Tumor','PanIN','Ductal_reactive','Duct_like_1','Duct_like_2','Acinar')
library(scCustomize)
polychrome_pal <- DiscretePalette_scCustomize(num_colors = length(unique(obj$cell_type_merge_v3_reduced)), palette = "polychrome")
coordinates <- Embeddings(obj, reduction = "umap.50PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length
p2 <- rasterize(DimPlot(obj, group.by = "cell_type_merge_v3_reduced", reduction = "umap.50PC", cols = polychrome_pal, pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + theme(legend.position="none") + ggtitle(NULL) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
p5 <- rasterize(DimPlot(obj, group.by = "cell_type_merge_v3_reduced", reduction = "umap.50PC", cols = polychrome_pal, pt.size=0.3, label=F, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
pdf(paste0(out_dir,"/",prefix,"_cell_type_merge_v3_reduced_fig6.pdf"), useDingbats = F,height=5, width = 5)
print(p2)
print(p5)
dev.off()

# extracting all cell types to a single table for the supplementary information.
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/
# conda activate seurat
library(tidyverse)
set.seed(1234)
compiling_input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/compiling_input_table_v7.tsv",sep = '\t',header=T)
sample_vector <- compiling_input_table$sample_ID
for (i in 1:length(sample_vector)) {
    sample_ID = sample_vector[i]
    print(sample_ID)
    is_merged = compiling_input_table[(compiling_input_table$sample_ID == sample_ID),"included_in_merged_object"]
    barcode_path <- compiling_input_table[(compiling_input_table$sample_ID == sample_ID),"Path_to_final_cell_type"]
    individual_path <- compiling_input_table[(compiling_input_table$sample_ID == sample_ID),"Path_to_original_individual_cell_type"]
    if (i == 1) {
        cell_barcode_table <- read.table(barcode_path, sep = '\t', header = T)
        colnames(cell_barcode_table) <- c("Xenium_barcode","cell_type_v3")
        cell_barcode_table$Sample_ID <- sample_ID
        cell_barcode_table$included_in_merged_object <- is_merged
        cell_barcode_table$merge_object_xenium_barcode <- NA
        cell_barcode_table$cell_type_merge_v3_reduced <- NA
        if (is_merged) {
            cell_barcode_table$merge_object_xenium_barcode[cell_barcode_table$Sample_ID == sample_ID] <- paste0(sample_ID,"_",cell_barcode_table$Xenium_barcode)
            cell_barcode_table$cell_type_merge_v3_reduced <- cell_barcode_table$cell_type_v3
            cell_barcode_table$cell_type_merge_v3_reduced[cell_barcode_table$cell_type_merge_v3 %in% c('CD4+T_cell','CD8+T_cell','Treg','Naive_T_cell')] <- "T_cell"
            cell_barcode_table$cell_type_merge_v3_reduced[cell_barcode_table$cell_type_merge_v3 %in% c('myCAF','iCAF')] <- "Fibroblast"
            cell_barcode_table$cell_type_merge_v3_reduced[cell_barcode_table$cell_type_merge_v3 %in% c('Islet_alpha','Islet_beta','Islet_delta','Islet_gamma','Islet_INS+GCG+')] <- "Islet"
            cell_barcode_table$cell_type_merge_v3_reduced[cell_barcode_table$cell_type_merge_v3 %in% c('Ductal_reactive','Duct_like_1','Duct_like_2')] <- "Ductal"
        }
    } else {
        cell_barcodes <- read.table(barcode_path, sep = '\t', header = T)
        colnames(cell_barcodes) <- c("Xenium_barcode","cell_type_v3")
        cell_barcodes$Sample_ID <- sample_ID
        cell_barcodes$included_in_merged_object <- is_merged
        cell_barcodes$merge_object_xenium_barcode <- NA
        cell_barcodes$cell_type_merge_v3_reduced <- NA
        if (is_merged) {
            cell_barcodes$merge_object_xenium_barcode[cell_barcodes$Sample_ID == sample_ID] <- paste0(sample_ID,"_",cell_barcodes$Xenium_barcode)
            cell_barcodes$cell_type_merge_v3_reduced <- cell_barcodes$cell_type_v3
            cell_barcodes$cell_type_merge_v3_reduced[cell_barcodes$cell_type_merge_v3 %in% c('CD4+T_cell','CD8+T_cell','Treg','Naive_T_cell')] <- "T_cell"
            cell_barcodes$cell_type_merge_v3_reduced[cell_barcodes$cell_type_merge_v3 %in% c('myCAF','iCAF')] <- "Fibroblast"
            cell_barcodes$cell_type_merge_v3_reduced[cell_barcodes$cell_type_merge_v3 %in% c('Islet_alpha','Islet_beta','Islet_delta','Islet_gamma','Islet_INS+GCG+')] <- "Islet"
            cell_barcodes$cell_type_merge_v3_reduced[cell_barcodes$cell_type_merge_v3 %in% c('Ductal_reactive','Duct_like_1','Duct_like_2')] <- "Ductal"
        }
        print(head(cell_barcodes))
        cell_barcode_table <- rbind(cell_barcode_table,cell_barcodes)
    }
}
write.table(cell_barcode_table, "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/all_cell_types_v7.tsv",sep='\t',quote=F)
for (i in 1:length(sample_vector)) {
    sample_ID = sample_vector[i]
    print(sample_ID)
    sample_barcode_table <- cell_barcode_table[(cell_barcode_table$Sample_ID == sample_ID), c("Xenium_barcode","cell_type_v3")]
    print(dim(sample_barcode_table))
    colnames(sample_barcode_table) <- c("cell_id","group")
    write.table(sample_barcode_table, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/",sample_ID,"_cell_types_v7.tsv"),sep='\t',quote=F,row.names=F)
    write.table(sample_barcode_table, paste0("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/",sample_ID,"_cell_types_v7.csv"),sep=',',quote=F,row.names=F)
}

#!/usr/bin/env Rscript
# Plotting the Xenium QC FeaturePlots
# supplementary figure control probes
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(ggh4x)
options(future.globals.maxSize= +Inf) # need to increase for SP001P1 #+Inf is 100% overkill but I don't care.
set.seed(1234)
# all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
# input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
# out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table <- read.table(args[1] ,sep='\t',header=T)
out_dir = args[2]
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
    DefaultAssay(obj) <- "Xenium.with.snvs" # this needs to be done befor you crop or the crop will not be associated with the correct assay and you can't plot variant probes
    Idents(obj) <- "seurat_clusters"
    DefaultFOV(obj, assay='Xenium.with.snvs') <- 'fov.with.snvs'
    DefaultBoundary(obj[["fov.with.snvs"]]) <- "segmentation"
    cell_types <- read.table(Manual_cell_type, header = T, sep = '\t')
    obj$barcode <- rownames(obj@meta.data)
    colnames(cell_types) <- c("barcode","cell_type")
    rownames(cell_types) <- cell_types$barcode
    cell_types$barcode <- NULL
    obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
    # since the cell types were generated on an object that was subset there are some 
    # cells that were removed when we filtered certain genes out of the matrix prior to integration.
    obj$cell_type[is.na(obj$cell_type)] <- "LowCount"
    unknown_low_quality_labels_v3 <- unique(unknown_low_quality_labels_v3,"LowCount")
    obj$neoplasm_normal_unknown <- NA
    obj$neoplasm_normal_unknown[(obj$cell_type  %in% tumor_labels)] <- "Neoplastic"
    #obj$neoplasm_normal_unknown[(obj$cell_type == "PanIN" )] <- "PanIN"
    obj$neoplasm_normal_unknown[(!(obj$cell_type %in% tumor_labels))] <- "Normal"
    #obj$neoplasm_normal_unknown[((obj$cell_type %in% c("Duct_like_1","Duct_like_2")))] <- "normal_duct"
    obj$neoplasm_normal_unknown[(obj$cell_type %in% unknown_low_quality_labels_v3)] <- "LowCount/Unknown"
    #obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("normal","normal_duct","PDAC","PanIN","low_quality/unknown"))
    obj$neoplasm_normal_unknown <- factor(obj$neoplasm_normal_unknown, levels = c("Normal","Neoplastic","LowCount/Unknown"))
    coordinates <- Embeddings(obj, reduction = "umap.30PC")
    dim_names <- colnames(coordinates)
    dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
    dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
    sorted_dim_ylims <- sort(dim_ylims)
    sorted_dim_xlims <- sort(dim_xlims)
    x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
    y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
    xy_ratio <- x_length/y_length
    print(paste0("plotting_",sample))
    pdf(paste0(out_dir,"/",sample,"/FeaturePlot_QC_",sample,".pdf"),useDingbats = F, height=6, width = 7)
    p1 <- FeaturePlot(obj, features = "nCount_BlankCodeword", order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p2 <- FeaturePlot(obj, features = "nFeature_BlankCodeword", order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p3 <- FeaturePlot(obj, features = "nCount_ControlCodeword", order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p4 <- FeaturePlot(obj, features = "nFeature_ControlCodeword", order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p5 <- FeaturePlot(obj, features = "nCount_ControlProbe", order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p6 <- FeaturePlot(obj, features = "nFeature_ControlProbe", order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p7 <- FeaturePlot(obj, features = "nCount_ControlProbe", max.cutoff = 1, order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    p8 <- FeaturePlot(obj, features = "nFeature_ControlProbe", max.cutoff = 1, order = T, reduction = "umap.30PC", label=F, pt.size=0.2, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio) + force_panelsizes(rows = unit(4, "in"), cols = unit(4, "in"))
    print(rasterize(p1, layers='Point', dpi=300))
    print(rasterize(p2, layers='Point', dpi=300))
    print(rasterize(p3, layers='Point', dpi=300))
    print(rasterize(p4, layers='Point', dpi=300))
    print(rasterize(p5, layers='Point', dpi=300))
    print(rasterize(p6, layers='Point', dpi=300))
    p7 <- ImageFeaturePlot(obj, features = "nCount_BlankCodeword", size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p8 <- ImageFeaturePlot(obj, features = "nFeature_BlankCodeword", size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p9 <- ImageFeaturePlot(obj, features = "nCount_ControlCodeword", size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p10 <- ImageFeaturePlot(obj, features = "nFeature_ControlCodeword", size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p11 <- ImageFeaturePlot(obj, features = "nCount_ControlProbe", size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p12 <- ImageFeaturePlot(obj, features = "nFeature_ControlProbe", size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p13 <- ImageFeaturePlot(obj, features = "nCount_ControlProbe",  max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p14 <- ImageFeaturePlot(obj, features = "nFeature_ControlProbe",  max.cutoff = 1, size = 0.15, cols = c("darkblue", "yellow"), border.size= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    pdf(paste0(out_dir,"/",sample,"/ImageFeaturePlot_QC_Control_Probes_",sample,".pdf"),useDingbats = F, height=6, width = 7)
    print(rasterize(p7, layers='Polygon', dpi=600))
    print(rasterize(p8, layers='Polygon', dpi=600))
    print(rasterize(p9, layers='Polygon', dpi=600))
    print(rasterize(p10, layers='Polygon', dpi=600))
    print(rasterize(p11, layers='Polygon', dpi=600))
    print(rasterize(p12, layers='Polygon', dpi=600))
    print(rasterize(p13, layers='Polygon', dpi=600))
    print(rasterize(p14, layers='Polygon', dpi=600))
    dev.off()
}


