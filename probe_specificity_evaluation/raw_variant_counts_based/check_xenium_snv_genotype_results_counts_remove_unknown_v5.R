# counts
library(Seurat)
library(tidyverse)
library(optparse)
library(epitools)
library(effectsize)
set.seed(1234)
# Rscript /diskmnt/Projects/Users/austins2/tools/xenium_snvs/check_xenium_snv_genotype_results_counts_remove_unknown_v5.R \
# -i /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/individual_objects/HT242P1-S1H4L4U1/HT242P1-S1H4L4U1_processed.rds \
# -o /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6/counts_based/HT242P1-S1H4L4U1/ \
# -c /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/compiled_cell_types/HT242P1-S1H4L4U1_cell_types_v7.tsv \
# -v /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv \
# -r /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_reference_allele_probes_hyphens2.tsv \
# -s HT242P1-S1H4L4U1 \
# -t 'Tumor' \
# -u 'LowCount'
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
    make_option(c("-u", "--unknown_cell_label"),
                type="character",
                default="",
                help="comma separated list of unknown/necrotic cell labels in the object. All cells with these labels will be excluded completely from the contingency table and significance testing calculations.",
                metavar="character")
);
# read the inputs and create the outdir
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
obj <- readRDS(opt$input)
cell_types <- read.table(opt$cell_types, header = T, sep = '\t')
out_dir <- paste0(opt$output,"/")
dir.create(out_dir, showWarnings = F, recursive = T)
tumor_cell_labels <- opt$tumor_cell_label
tumor_cell_labels <- unlist(strsplit(tumor_cell_labels, ","))
unknown_labels <- opt$unknown_cell_label
unknown_labels <- unlist(strsplit(unknown_labels, ","))
setwd(out_dir)
variant_alleles <- read.table(opt$variant_alleles,sep='\t',header=F)
#print(variant_alleles)
reference_alleles <- read.table(opt$reference_alleles,sep='\t',header=F)
out_sink <- file(paste0(out_dir,opt$sample,"_genotyping_results.txt"), open = "wt")
sink(file = out_sink, append = FALSE, type = "output", split = TRUE)
# sink(file = out_sink, append = FALSE, type = "message", split = FALSE)
# add the cell type labels to the object metadata
colnames(cell_types) <- c("barcode","cell_type")
rownames(cell_types) <- cell_types$barcode
cell_types$barcode <- NULL
obj <- AddMetaData(object = obj, metadata = cell_types, col.name = "cell_type")
#print("removing the following cell types before calculation:")
#print(unknown_labels)
cell_type_vector <- unique(cell_types$cell_type)
retain_cell_types <- cell_type_vector[(!(cell_type_vector %in% unknown_labels))]
#print(retain_cell_types)
obj <- subset(obj, subset = cell_type %in% retain_cell_types)
#print("proceeding with just the following cell types:")
#print(unique(obj$cell_type))
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
for (assay in assays_with_snv_probes) {
    number_of_variants <- length(variant_alleles[,"V1"])
    statistics_table_results <- data.frame(Variant.probe.name = rep(NA, number_of_variants),
                                           Reference.probe.name = rep(NA, number_of_variants),
                                           Number.of.variant.probes.in.cancer.cells = rep(NA, number_of_variants),
                                           Number.of.reference.probes.in.cancer.cells = rep(NA, number_of_variants),
                                           Number.of.variant.probes.in.non.cancer.cells = rep(NA, number_of_variants),
                                           Number.of.reference.probes.in.non.cancer.cells = rep(NA, number_of_variants),
                                           total.number.of.both.variant.site.probes.in.cancer.cells = rep(NA, number_of_variants),
                                           total.number.of.both.variant.site.probes.in.normal.cells = rep(NA, number_of_variants),
                                           proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = rep(NA, number_of_variants), # Proportion.test.value.1 # true positive rate
                                           proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = rep(NA, number_of_variants), # Proportion.test.value.2
                                           Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = rep(NA, number_of_variants),
                                           Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = rep(NA, number_of_variants),
                                           Chi.square.p.value = rep(NA, number_of_variants),
                                           Proportion.test.95CI.low = rep(NA, number_of_variants),
                                           Proportion.test.95CI.high = rep(NA, number_of_variants),
                                           Proportion.test.p.value = rep(NA, number_of_variants),
                                           fisher.oddsratio = rep(NA, number_of_variants),
                                           fisher.oddsratio.95CI.low = rep(NA, number_of_variants),
                                           fisher.oddsratio.95CI.high = rep(NA, number_of_variants),
                                           pearsons.phi.effectsize = rep(NA, number_of_variants),
                                           pearsons.phi.effectsize.95CI.low = rep(NA, number_of_variants),
                                           pearsons.phi.effectsize.95CI.high = rep(NA, number_of_variants),
                                           bonferroni.corrected.chisq.test.pval = rep(NA, number_of_variants),
                                           bonferroni.corrected.prop.test.pval = rep(NA, number_of_variants),
                                           total.number.of.cancer.cells.in.sample = rep(NA, number_of_variants),
                                           total.number.of.normal.cells.in.sample = rep(NA, number_of_variants),
                                           total.number.of.cells.in.sample.after.low.quality.removed = rep(NA, number_of_variants),
                                           row.names = c(1:number_of_variants))
    print(paste0("SNV probes found in ",assay," assay. Currently testing ",assay," assay."))
    assay_class = class(obj[[assay]])
    if (assay_class == "Assay5") {
        #print(1)
        counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
        counts_df <- as.data.frame(counts_df)
        #print(counts_df[1:5,1:5])
    } else {
        counts_df <- as.data.frame(t(as.matrix(GetAssayData(object = obj, assay = assay, slot = "counts"))))
    }
    for (i in 1:length(variant_alleles[,"V1"])) {
        variant_key <- variant_alleles[i,'V1']
        wt_key <- reference_alleles[i,'V1']
        statistics_table_results[i, "Variant.probe.name"] <- variant_key
        statistics_table_results[i, "Reference.probe.name"] <- wt_key
        print(variant_key)
        print(wt_key)
        alt_tmp <- as.data.frame(counts_df[ ,variant_key],row.names=rownames(counts_df))
        wt_tmp <- as.data.frame(counts_df[ ,wt_key], row.names=rownames(counts_df))
        obj <- AddMetaData(object = obj, metadata = alt_tmp, col.name = paste0(variant_key,"_",assay,"_count"))
        obj <- AddMetaData(object = obj, metadata = wt_tmp, col.name = paste0(wt_key,"_",assay,"_count"))
        obj@meta.data[paste0(variant_key,"_",assay,"_call")] <- NA
        obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(wt_key,"_",assay,"_count")] > 0] <- "WT"
        obj@meta.data[paste0(variant_key,"_",assay,"_call")][obj@meta.data[paste0(variant_key,"_",assay,"_count")] > 0] <- "ALT"
        obj$tumor_cell_call <- "Normal_cell" # compress the cell type labels to just normal or tumor
        #obj$tumor_cell_call[obj$cell_type %in% c("Tumor","Cancer","Cancer cells","Cancer cells proliferating","PC","ITPN","Tumor_proliferating","Tumor_Proliferative")] <- "Cancer_cell" # compress the cell type labels to just normal or tumor
        obj$tumor_cell_call[obj$cell_type %in% tumor_cell_labels] <- "Cancer_cell" # compress the cell type labels to just normal or tumor
        statistics_table_results[i, "total.number.of.cancer.cells.in.sample"] <- sum(obj$tumor_cell_call == "Cancer_cell")
        statistics_table_results[i, "total.number.of.normal.cells.in.sample"] <- sum(obj$tumor_cell_call == "Normal_cell")
        statistics_table_results[i, "total.number.of.cells.in.sample.after.low.quality.removed"] <- sum(obj$tumor_cell_call == "Normal_cell") + sum(obj$tumor_cell_call == "Cancer_cell")
        #data <- table(obj@meta.data$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")])
        data <- data.frame( ALT = c(0,0), # Cancer_cell, Normal_Cell
                            WT = c(0,0), # Cancer_cell, Normal_Cell
                            row.names = c("Cancer_cell", "Normal_cell"))
        data["Normal_cell","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"),paste0(wt_key,"_",assay,"_count")])
        data["Normal_cell","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal_cell"),paste0(variant_key,"_",assay,"_count")])
        data["Cancer_cell","WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"),paste0(wt_key,"_",assay,"_count")])
        data["Cancer_cell","ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_cell"),paste0(variant_key,"_",assay,"_count")])
        statistics_table_results[i, "Number.of.variant.probes.in.cancer.cells"] <- data["Cancer_cell","ALT"]
        statistics_table_results[i, "Number.of.reference.probes.in.cancer.cells"] <- data["Cancer_cell","WT"]
        statistics_table_results[i, "Number.of.variant.probes.in.non.cancer.cells"] <- data["Normal_cell","ALT"]
        statistics_table_results[i, "Number.of.reference.probes.in.non.cancer.cells"] <- data["Normal_cell","WT"]
        statistics_table_results[i, "total.number.of.both.variant.site.probes.in.cancer.cells"] <- data["Cancer_cell","ALT"] + data["Cancer_cell","WT"]
        statistics_table_results[i, "total.number.of.both.variant.site.probes.in.normal.cells"] <- data["Normal_cell","ALT"] + data["Normal_cell","WT"]
        if ((sum(counts_df[ ,variant_key]) > 0) & (sum(counts_df[ ,wt_key]) >0)) {
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
                    oddsratio.out = oddsratio(as.matrix(data), method = "fisher")
                    print("Odds ratio")
                    print(oddsratio.out)
                    phi.effectsize <- phi(as.matrix(data), digits = 3)
                    print("Effect size (Pearson's phi)")
                    print(phi.effectsize)
                    #variant_fisherexact_test <- fisher.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE,alternative="less")
                    #print("simulate.p.value=FALSE", "alternative=less")
                    #print(variant_fisherexact_test)
                    #variant_fisherexact_test <- fisher.test(obj$tumor_cell_call, obj@meta.data[,paste0(variant_key,"_",assay,"_call")], simulate.p.value=FALSE,alternative="two.sided")
                    #print("simulate.p.value=FALSE", "alternative=two.sided")
                    #print(variant_fisherexact_test)
                    
                    statistics_table_results[i, 'proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells'] <- variant_prop_test$estimate[1][[1]] # true.positive.rate if the sample actually has the mutation
                    statistics_table_results[i, 'proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells'] <- variant_prop_test$estimate[2][[1]] # false.positive.rate if the sample actually has the mutation.
                    statistics_table_results[i, 'Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells'] <- data["Cancer_cell","WT"] / (data["Cancer_cell","ALT"] + data["Cancer_cell","WT"]) # false.negative.rate if the sample actually has the mutation
                    statistics_table_results[i, 'Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells'] <- data["Normal_cell","WT"] / (data["Normal_cell","ALT"] + data["Normal_cell","WT"]) # true.negative.rate if the sample actually has the mutation
                    
                    statistics_table_results[i, 'Chi.square.p.value'] <- variant_chisq_test1$p.value
                    statistics_table_results[i, 'bonferroni.corrected.chisq.test.pval'] <- p.adjust(variant_chisq_test1$p.value, method = "bonferroni", n = number_of_variants)
                    statistics_table_results[i, 'Proportion.test.p.value'] <- variant_prop_test$p.value
                    statistics_table_results[i, 'bonferroni.corrected.prop.test.pval'] <- p.adjust(variant_prop_test$p.value, method = "bonferroni", n = number_of_variants)
                    statistics_table_results[i, 'Proportion.test.95CI.low'] <- variant_prop_test$conf.int[1]
                    statistics_table_results[i, 'Proportion.test.95CI.high'] <- variant_prop_test$conf.int[2]
                    statistics_table_results[i, 'fisher.oddsratio'] <- oddsratio.out$Odds_ratio
                    statistics_table_results[i, 'fisher.oddsratio.95CI.low'] <- oddsratio.out$CI_low
                    statistics_table_results[i, 'fisher.oddsratio.95CI.high'] <- oddsratio.out$CI_high
                    statistics_table_results[i, 'pearsons.phi.effectsize'] <- phi.effectsize$phi_adjusted
                    statistics_table_results[i, 'pearsons.phi.effectsize.95CI.low'] <- phi.effectsize$CI_low
                    statistics_table_results[i, 'pearsons.phi.effectsize.95CI.high'] <- phi.effectsize$CI_high
                } else {
                    print(paste0(sum(data[ "Cancer_cell",])," counts of the variant and wildtype allele were detected unable to run a chi-square or proportion test"))
                    print(paste0(sum(data[ "Normal_cell",])," counts of the variant and wildtype allele were detected unable to run a chi-square or proportion test"))
                    print(paste0(sum(data["Cancer_cell","ALT"])," counts of the variant allele in Cancer cells were detected unable to run a chi-square or proportion test"))
                    print(paste0(sum(data["Normal_cell","ALT"])," counts of the variant allele Normal cells were detected unable to run a chi-square or proportion test"))
                    print(paste0(sum(data["Cancer_cell","WT"])," counts of the wildtype allele in Cancer cells were detected unable to run a chi-square or proportion test"))
                    print(paste0(sum(data["Normal_cell","WT"])," counts of the wildtype allele Normal cells were detected unable to run a chi-square or proportion test"))
                }
            } else {
                if ("Normal_cell" %in% rownames(data)) {
                    print("No counts of either allele in cancer cells ",sum(data["Normal_cell","ALT"])," counts of the variant allele in Normal cells ",sum(data["Normal_cell","WT"])," counts of the wildtype allele in Normal cells were detected unable to run a chi-square or proportion test")
                } else {
                    print("No counts of either allele in Normal cells ",sum(data["Cancer_cell","ALT"])," counts of the variant allele in Normal cells ",sum(data["Cancer_cell","WT"])," counts of the wildtype allele in Cancer cells were detected unable to run a chi-square or proportion test")
                }
            }
        } else {
            print(paste0(sum(counts_df[ ,variant_key])," counts of the variant allele were detected unable to run a chi-square or proportion test"))
            print(paste0(sum(counts_df[ ,wt_key])," counts of the wildtype allele were detected unable to run a chi-square or proportion test"))
        }
    }
    #statistics_table_results$fisher.oddsratio[is.infinite(statistics_table_results$fisher.oddsratio)] = NA
    #statistics_table_results$pearsons.phi.effectsize[is.infinite(statistics_table_results$pearsons.phi.effectsize)] = NA
    write.table(statistics_table_results, paste0(out_dir, opt$sample, "_probe_specificity_table.csv"),quote=F,sep=',')
}
sink()
# sink()
close(out_sink)
write.table(obj@meta.data,paste0(out_dir,opt$sample,"_genotyping_metadata.tsv"),quote=F,sep='\t')
