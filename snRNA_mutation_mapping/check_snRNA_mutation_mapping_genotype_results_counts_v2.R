library(Seurat)
library(tidyverse)
library(optparse)
library(ggrastr)
set.seed(1234)
# Example command:
# Rscript /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/mutation_mapping/snRNA_mapping_freeze_v1_2024-05-22/proportion_tests/cell_type_v1/check_snRNA_mutation_mapping_genotype_results_counts_v2.R -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1_processed_update_no_doublet_seurat5.0.1.rds -o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/mutation_mapping/snRNA_mapping_freeze_v3_2024-11-25_c3.1/proportion_tests/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1 -c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/matching_snRNA_objects/seurat_v5.0.1/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1_cell_type_xenium_variant_v1.tsv -v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv -r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv -m /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/mutation_mapping/matrix_maf_to_probe_feature_conversion_table.tsv -s HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1 -t 'Tumor,PanIN,Tumor_proliferating' -w /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/mutation_mapping/snRNA_mapping_freeze_v3_2024-11-25/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1.parsed.reference_allele_count_matrix.tsv -a /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/mutation_mapping/snRNA_mapping_freeze_v3_2024-11-25/HT227P1-S1H1A3Y1Nd1_1Z1_1Bmn1_1.parsed.variant_allele_count_matrix.tsv
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL, 
                help="path to input rds objects",
                metavar="character"),
    make_option(c("-o", "--output"),
                type="character",
                default="./", 
                help="output folder path",
                metavar="character"),
    make_option(c("-c", "--cell_types"),
                type="character",
                default="", 
                help="path to tsv file with cell types for the object. All cells must have a cell type. the cell barcode needs to be the first column. the cell type needs to be the second column. The columns need to have a header.",
                metavar="character"),
    make_option(c("-v", "--variant_alleles"),
                type="character",
                default="", 
                help="tsv file with the variant probes in the first column",
                metavar="character"),
    make_option(c("-r", "--reference_alleles"),
                type="character",
                default="", 
                help="tsv file with the reference probes in the first column",
                metavar="character"),
    make_option(c("-s", "--sample"),
                type="character",
                default="test", 
                help="name of the sample for the test. Output file names will contain the sample.",
                metavar="character"),
    make_option(c("-t", "--tumor_cell_label"),
                type="character",
                default="Tumor,Cancer,Cancer cells,Cancer cells proliferating,PC,ITPN,Tumor_proliferating,Tumor_Proliferative",
                help="Tumor cell labels in the object. All cells with these labels will be included in the 'Cancer_cell' label in the contingency table.",
                metavar="character"),
    make_option(c("-m", "--matrix_labels"),
                type="character",
                default="", 
                help="tsv file with the matrix labels in the first column, the variant probes in the second column, the reference probes in the third column",
                metavar="character"),
    make_option(c("-w", "--wildtype_allele_counts"),
                type="character",
                default="", 
                help="tsv matrix with counts of the wildtype allele in each cell",
                metavar="character"),
    make_option(c("-a", "--alternate_allele_counts"),
                type="character",
                default="", 
                help="tsv matrix with counts of the alternate allele in each cell",
                metavar="character")
);
# read the inputs and create the outdir
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
obj <- readRDS(opt$input)
cell_types <- read.table(opt$cell_types, header = F, sep = '\t')
out_dir <- paste0(opt$output,"/")
dir.create(out_dir, showWarnings = F, recursive = T)
tumor_cell_labels <- opt$tumor_cell_label
tumor_cell_labels <- unlist(strsplit(tumor_cell_labels, ","))
setwd(out_dir)
sample <- opt$sample
variant_alleles <- read.table(opt$variant_alleles,sep='\t',header=F)
#print(variant_alleles)
reference_alleles <- read.table(opt$reference_alleles,sep='\t',header=F)
#####
matrix_labels <- read.table(opt$matrix_labels,sep='\t', header=T, row.names="X")
reference_counts <- read.table(opt$wildtype_allele_counts,sep='\t', header=T, row.names="X0")
variant_counts <- read.table(opt$alternate_allele_counts,sep='\t', header=T, row.names="X0")
#####
out_sink <- file(paste0(out_dir,opt$sample,"_genotyping_results.txt"), open = "wt")
sink(file = out_sink, append = FALSE, type = "output", split = TRUE)
# sink(file = out_sink, append = FALSE, type = "message", split = FALSE)
# add the cell type labels to the object metadata
colnames(cell_types) <- c("barcode","cell_type")
cell_types$original_RNA_barcode <- stringr::str_sub(cell_types$barcode, -18, -1)
#####
obj$barcode <- rownames(obj@meta.data)
obj$original_RNA_barcode <- stringr::str_sub(rownames(obj@meta.data), -18, -1)
barcode_1 <- obj$barcode[1][[1]]
original_RNA_barcode_1 <- obj$original_RNA_barcode[1][[1]]
barcode_prefix = ""
if (!(barcode_1 == original_RNA_barcode_1)) {
    barcode_prefix = stringr::str_sub(barcode_1, 1, nchar(barcode_1)-18)
}
rownames(cell_types) <- paste0(barcode_prefix,cell_types$original_RNA_barcode)
cell_types$barcode <- NULL
cell_types$original_RNA_barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
#####
assay_list <- Assays(obj)
assays_with_snv_probes <- c()
#print("variant_alleles[1,'V1']")
#print(variant_alleles[1,'V1'])
for (assay in assay_list) {
    DefaultAssay(obj) <- assay
    feature_vector <- dimnames(obj)[1][[1]]
    #print("feature_vector")
    #print(feature_vector)
    if (variant_alleles[1,'V1'] %in% feature_vector) {
        assays_with_snv_probes <- c(assays_with_snv_probes, assay)
    }
}
print("assays with snv probes:")
print(assays_with_snv_probes)
# for (assay in assays_with_snv_probes) {
#     print(paste0("SNV probes found in ",assay," assay. Currently testing ",assay," assay."))
#     assay_class = class(obj[[assay]])
#     if (assay_class == "Assay5") {
#         #print(1)
#         counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
#         counts_df <- as.data.frame(counts_df)
#         #print(counts_df[1:5,1:5])
#     } else {
#         counts_df <- as.data.frame(t(as.matrix(GetAssayData(object = obj, assay = assay, slot = "counts"))))
#     }
number_of_variants = length(matrix_labels[,"matrix_key"])
contingency_table_results <- data.frame(Alternate.probe.name = matrix_labels[,"alt_probe_feature"], 
                                        Wild.type.probe.name = matrix_labels[,'ref_probe_feature'],
                                        snRNA_mutation_mapping_key = matrix_labels[,"matrix_key"],
                                        Chi.square.p.value = rep(NA, number_of_variants),
                                        Proportion.test.p.value = rep(NA, number_of_variants),
                                        Proportion.test.lower.interval.value = rep(NA, number_of_variants),
                                        Proportion.test.upper.interval.value = rep(NA, number_of_variants),
                                        Proportion.test.value.1 = rep(NA, number_of_variants),
                                        Proportion.test.value.2 = rep(NA, number_of_variants),
                                        Number.of.cancer.cells.with.the.alternate.probe = rep(NA, number_of_variants),
                                        Number.of.cancer.cells.with.the.wild.type.probe = rep(NA, number_of_variants),
                                        Number.of.non.cancer.cells.with.the.alternate.probe = rep(NA, number_of_variants),
                                        Number.of.non.cancer.cells.with.the.wild.type.probe = rep(NA, number_of_variants),
                                        sample.ID = rep(sample, number_of_variants),
                                        bonferroni.corrected.chisq.test.significance = rep(NA, number_of_variants),
                                        bonferroni.corrected.prop.test.significance = rep(NA, number_of_variants),
                                        row.names = c(1:number_of_variants))
print("SNV probes not present in assay. Adding mutation mapping matrix with counts of varaiant and wt alleles.")
pdf(paste0(sample,"_genotype_featureplots.pdf"), useDingbats=F, height=5, width = 7)
for (i in 1:number_of_variants) {
    variant_key <- matrix_labels[i,'alt_probe_feature']
    wt_key <- matrix_labels[i,'ref_probe_feature']
    #####
    matrix_key <- matrix_labels[i,"matrix_key"]
    #####
    print(variant_key)
    print(wt_key)
    #alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
    #wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
    wt_counts_df <- as.data.frame(reference_counts[,matrix_key], row.names = paste0(barcode_prefix, rownames(reference_counts)))
    alt_counts_df <- as.data.frame(variant_counts[,matrix_key], row.names = paste0(barcode_prefix, rownames(variant_counts)))
    colnames(wt_counts_df) <- c(matrix_key)
    colnames(alt_counts_df) <- matrix_key
    #if ((sum(counts_df[ ,variant_key]) > 0) & (sum(counts_df[ ,wt_key]) >0)) {
    obj <- AddMetaData(object = obj, metadata = alt_counts_df, col.name = paste0(matrix_key,"_",assay,"_alternate_count"))
    obj <- AddMetaData(object = obj, metadata = wt_counts_df, col.name = paste0(matrix_key,"_",assay,"_wildtype_count"))
    obj@meta.data[paste0(matrix_key,"_",assay,"_call")] <- NA
    obj@meta.data[paste0(matrix_key,"_",assay,"_call")][obj@meta.data[paste0(matrix_key,"_",assay,"_wildtype_count")] > 0] <- "WT"
    obj@meta.data[paste0(matrix_key,"_",assay,"_call")][obj@meta.data[paste0(matrix_key,"_",assay,"_alternate_count")] > 0] <- "ALT"
    obj$tumor_cell_call <- "Normal_cell" # compress the cell type labels to just normal or tumor
    obj$tumor_cell_call[obj$cell_type %in% tumor_cell_labels] <- "Cancer_cell"
    data <- as.data.frame(table(obj@meta.data$tumor_cell_call, obj@meta.data[,paste0(matrix_key,"_",assay,"_call")]))
    data <- data.frame( ALT = c(sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"), paste0(matrix_key,"_",assay,"_alternate_count")]),
                                sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"), paste0(matrix_key,"_",assay,"_alternate_count")])),
                        WT = c(sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"), paste0(matrix_key,"_",assay,"_wildtype_count")]),
                               sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"), paste0(matrix_key,"_",assay,"_wildtype_count")])),
                        row.names = c("Cancer_cell", "Normal_cell")
    )
    print(data)
    contingency_table_results[i, 'Number.of.non.cancer.cells.with.the.wild.type.probe'] <- data["Normal_cell","WT"]
    contingency_table_results[i, "Number.of.non.cancer.cells.with.the.alternate.probe"] <- data["Normal_cell","ALT"]
    contingency_table_results[i, "Number.of.cancer.cells.with.the.wild.type.probe"] <- data["Cancer_cell","WT"]
    contingency_table_results[i, "Number.of.cancer.cells.with.the.alternate.probe"] <- data["Cancer_cell","ALT"]
    print(paste0("generating plot of ",matrix_key))
    p1 <- FeaturePlot(obj, features = paste0(matrix_key,"_",assay,"_alternate_count"), order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + labs(title = paste0(variant_key,"_alternate_count"), subtitle = paste0(matrix_key,"_alternate_count")) + theme(plot.title = element_text(size=12, hjust = 0.5))
    print(rasterize(p1, layers='Point', dpi=300))
    p2 <- FeaturePlot(obj, features = paste0(matrix_key,"_",assay,"_wildtype_count"), order = T, reduction = "umap.30PC", label=F, pt.size=0.3, raster=FALSE) + labs(title = paste0(wt_key,"_wildtype_count"), subtitle = paste0(matrix_key,"_wildtype_count")) + theme(plot.title = element_text(size=12, hjust = 0.5))
    print(rasterize(p2, layers='Point', dpi=300))
    if ((sum(alt_counts_df[ ,matrix_key]) > 0) & (sum(wt_counts_df[ ,matrix_key]) > 0)) {
        # #print(head(alt_tmp))
        # #print(head(wt_tmp))
        # # obj <- AddMetaData(object = obj, metadata = alt_counts_df, col.name = paste0(variant_key,"_",assay,"_count"))
        # # obj <- AddMetaData(object = obj, metadata = wt_counts_df, col.name = paste0(wt_key,"_",assay,"_count"))
        # obj <- AddMetaData(object = obj, metadata = alt_counts_df, col.name = paste0(matrix_key,"_",assay,"_alternate_count"))
        # obj <- AddMetaData(object = obj, metadata = wt_counts_df, col.name = paste0(matrix_key,"_",assay,"_wildtype_count"))
        # # obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
        # obj@meta.data[paste0(matrix_key,"_",assay,"_call")] <- NA
        # # obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
        # # obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
        # obj@meta.data[paste0(matrix_key,"_",assay,"_call")][obj@meta.data[paste0(matrix_key,"_",assay,"_wildtype_count")] > 0] <- "WT"
        # obj@meta.data[paste0(matrix_key,"_",assay,"_call")][obj@meta.data[paste0(matrix_key,"_",assay,"_alternate_count")] > 0] <- "ALT"
        # obj$tumor_cell_call <- "Normal_cell" # compress the cell type labels to just normal or tumor
        # #obj$tumor_cell_call[obj$cell_type %in% c("Tumor","Cancer","Cancer cells","Cancer cells proliferating","PC","ITPN","Tumor_proliferating","Tumor_Proliferative")] <- "Cancer_cell" # compress the cell type labels to just normal or tumor
        # obj$tumor_cell_call[obj$cell_type %in% tumor_cell_labels] <- "Cancer_cell" # compress the cell type labels to just normal or tumor
        # data <- table(obj@meta.data$tumor_cell_call, obj@meta.data[,paste0(matrix_key,"_",assay,"_call")])
        # # data <- table(obj@meta.data$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")])
        # # data["Normal_cell","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"), paste0(wt_key,"_",assay,"_count")])
        # # data["Normal_cell","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"),paste0(variant_key,"_",assay,"_count")])
        # # data["Cancer_cell","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"),paste0(wt_key,"_",assay,"_count")])
        # # data["Cancer_cell","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"),paste0(variant_key,"_",assay,"_count")])
        # data["Normal_cell","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"), paste0(matrix_key,"_",assay,"_wildtype_count")])
        # data["Normal_cell","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"),paste0(matrix_key,"_",assay,"_alternate_count")])
        # data["Cancer_cell","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"),paste0(matrix_key,"_",assay,"_wildtype_count")])
        # data["Cancer_cell","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"),paste0(matrix_key,"_",assay,"_alternate_count")])
        if (("Normal_cell" %in% rownames(data)) & ("Cancer_cell" %in% rownames(data))) {
            if ((sum(data["Normal_cell",]) !=0) & (sum(data["Cancer_cell",]) !=0)) {
                print(data)
                variant_chisq_test1 <- chisq.test(data, simulate.p.value=TRUE, correct=FALSE)
                #variant_chisq_test2 <- chisq.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=TRUE, correct=FALSE)
                print("simulate.p.value=TRUE, correct=FALSE")
                print(variant_chisq_test1)
                #print("simulate.p.value=TRUE, correct=FALSE")
                #print(variant_chisq_test2)
                #variant_chisq_test <- chisq.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE, correct=FALSE)
                #print("simulate.p.value=FALSE, correct=FALSE")
                #print(variant_chisq_test)
                #variant_chisq_test <- chisq.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE, correct=TRUE)
                #print("simulate.p.value=FALSE, correct=TRUE")
                #print(variant_chisq_test)
                #variant_fisherexact_test <- fisher.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE)#,alternative="greater")
                #print("simulate.p.value=FALSE")#, "alternative=greater")
                #print(variant_fisherexact_test)
                x <- data[,'ALT']
                #print(x)
                n <- rowSums(data)
                #print(n)
                variant_prop_test <- prop.test(x, n , alternative = "greater", correct = FALSE) 
                print("prop_test alternative=greater")
                print(variant_prop_test)
                #variant_fisherexact_test <- fisher.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE,alternative="less")
                #print("simulate.p.value=FALSE", "alternative=less")
                #print(variant_fisherexact_test)
                #variant_fisherexact_test <- fisher.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE,alternative="two.sided")
                #print("simulate.p.value=FALSE", "alternative=two.sided")
                #print(variant_fisherexact_test)
                contingency_table_results[i, 'Chi.square.p.value'] <- variant_chisq_test1$p.value
                contingency_table_results[i, 'bonferroni.corrected.chisq.test.significance'] <- p.adjust(variant_chisq_test1$p.value, method = "bonferroni", n = number_of_variants)
                contingency_table_results[i, 'Proportion.test.p.value'] <- variant_prop_test$p.value
                contingency_table_results[i, 'bonferroni.corrected.prop.test.significance'] <- p.adjust(variant_prop_test$p.value, method = "bonferroni", n = number_of_variants)
                contingency_table_results[i, 'Proportion.test.value.1'] <- variant_prop_test$estimate[1][[1]]
                contingency_table_results[i, 'Proportion.test.value.2'] <- variant_prop_test$estimate[2][[1]]
                contingency_table_results[i, 'Proportion.test.lower.interval.value'] <- variant_prop_test$conf.int[1]
                contingency_table_results[i, 'Proportion.test.upper.interval.value'] <- variant_prop_test$conf.int[2]
            } else {
                print(paste0(sum(data[ "Cancer_cell",])," counts of the alternate and wildtype allele were detected unable to run a chi-square or proportion test"))
                print(paste0(sum(data[ "Normal_cell",])," counts of the alternate and wildtype allele were detected unable to run a chi-square or proportion test"))
                print(paste0(sum(data["Cancer_cell","ALT"])," counts of the alternate allele in Cancer cells were detected unable to run a chi-square or proportion test"))
                print(paste0(sum(data["Normal_cell","ALT"])," counts of the alternate allele Normal cells were detected unable to run a chi-square or proportion test"))
                print(paste0(sum(data["Cancer_cell","WT"])," counts of the wildtype allele in Cancer cells were detected unable to run a chi-square or proportion test"))
                print(paste0(sum(data["Normal_cell","WT"])," counts of the wildtype allele Normal cells were detected unable to run a chi-square or proportion test"))
            }
        } else {
            if ("Normal_cell" %in% rownames(data)) {
                print("No counts of either allele in cancer cells ",sum(data["Normal_cell","ALT"])," counts of the alternate allele in Normal cells ",sum(data["Normal_cell","WT"])," counts of the wildtype allele in Normal cells were detected unable to run a chi-square or proportion test")
            } else {
                print("No counts of either allele in Normal cells ",sum(data["Cancer_cell","ALT"])," counts of the alternate allele in Normal cells ",sum(data["Cancer_cell","WT"])," counts of the wildtype allele in Cancer cells were detected unable to run a chi-square or proportion test")
            }
        }
    } else {
        print(paste0(sum(alt_counts_df[ ,matrix_key])," counts of the alternate allele were detected. unable to run a chi-square or proportion test"))
        print(paste0(sum(wt_counts_df[ ,matrix_key])," counts of the wildtype allele were detected. unable to run a chi-square or proportion test"))
        # print(paste0(sum(counts_df[ ,variant_key])," counts of the alternate allele were detected unable to run a chi-square or proportion test"))
        # print(paste0(sum(counts_df[ ,wt_key])," counts of the wildtype allele were detected unable to run a chi-square or proportion test"))
    }
}
dev.off()
sink()
# sink()
close(out_sink)
write.table(obj@meta.data,paste0(out_dir,opt$sample,"_genotyping_metadata.tsv"),quote=F,sep='\t')
write.table(contingency_table_results, paste0(sample, "_mapped_mutation_specificity_results_by_sample.tsv"), sep='\t',quote=F)
