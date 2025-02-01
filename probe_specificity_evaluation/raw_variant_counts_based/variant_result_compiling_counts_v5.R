library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table <- read.table(args[1],sep='\t',header=T,row.names="X")
out_dir = paste0(args[2],"/")
#variant_list <- read.table("/diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_SNV_panel_v1_Z3DKDB_variant_allele_probes_hyphens2.tsv", sep='\t', header=F)
# this contains the list of all variants present in the Xenium human custom panel v2 and the Xenium snv panel v1
variant_list <- read.table("/diskmnt/Projects/Users/austins2/tools/xenium_snvs/All_xenium_variant_allele_probes_hyphens.tsv", sep='\t', header=F)
dir.create(paste0(out_dir,"/Variant_specific_results/"))
variant_list <- variant_list$V1
# input_table <- read.table("plotting_input_table.tsv",sep='\t',header=T,row.names="X")
# # 1 sample_ID probe_specificity_test
iterables <- rownames(input_table)
probe_table_list <- list()
sample_vector <- input_table$sample_ID
print("reading input")
for (sample in sample_vector) {
    print(sample)
    print(input_table$probe_specificity_test[input_table$sample_ID == sample])
    probe_table_list[[sample]] <- read.table(input_table$probe_specificity_test[input_table$sample_ID == sample], sep=',',header=T)
}
saveRDS(probe_table_list,"probe_table_list.rds")
merged_probe_table <- data.frame(Variant.probe.name = c(),#probe_table_list[[sample_ID]][ ,]),
                                 Reference.probe.name = c(),
                                 Number.of.variant.probes.in.cancer.cells = c(),
                                 Number.of.reference.probes.in.cancer.cells = c(),
                                 Number.of.variant.probes.in.non.cancer.cells = c(),
                                 Number.of.reference.probes.in.non.cancer.cells = c(),
                                 total.number.of.both.variant.site.probes.in.cancer.cells = c(),
                                 total.number.of.both.variant.site.probes.in.normal.cells = c(),
                                 proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = c(), # Proportion.test.value.1 # true positive rate
                                 proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = c(), # Proportion.test.value.2
                                 Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = c(),
                                 Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = c(),
                                 Chi.square.p.value = c(),
                                 Proportion.test.95CI.low = c(),
                                 Proportion.test.95CI.high = c(),
                                 Proportion.test.p.value = c(),
                                 fisher.oddsratio = c(),
                                 fisher.oddsratio.95CI.low = c(),
                                 fisher.oddsratio.95CI.high = c(),
                                 pearsons.phi.effectsize = c(),
                                 pearsons.phi.effectsize.95CI.low = c(),
                                 pearsons.phi.effectsize.95CI.high = c(),
                                 sample.ID = c(),
                                 bonferroni.corrected.chisq.test.pval = c(),
                                 bonferroni.corrected.prop.test.pval = c(),
                                 total.number.of.cancer.cells.in.sample = c(),
                                 total.number.of.normal.cells.in.sample = c(),
                                 total.number.of.cells.in.sample.after.low.quality.removed = c(),
                                 stringsAsFactors = FALSE)
# sampleID probe_specificity_test
rename = 0
for (variant in variant_list) {
    print(variant)
    new_table <- data.frame(Variant.probe.name = c(),#probe_table_list[[sample_ID]][ ,]),
                            Reference.probe.name = c(),
                            Number.of.variant.probes.in.cancer.cells = c(),
                            Number.of.reference.probes.in.cancer.cells = c(),
                            Number.of.variant.probes.in.non.cancer.cells = c(),
                            Number.of.reference.probes.in.non.cancer.cells = c(),
                            total.number.of.both.variant.site.probes.in.cancer.cells = c(),
                            total.number.of.both.variant.site.probes.in.normal.cells = c(),
                            proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = c(), # Proportion.test.value.1 # true positive rate
                            proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = c(), # Proportion.test.value.2
                            Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = c(),
                            Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = c(),
                            Chi.square.p.value = c(),
                            Proportion.test.95CI.low = c(),
                            Proportion.test.95CI.high = c(),
                            Proportion.test.p.value = c(),
                            fisher.oddsratio = c(),
                            fisher.oddsratio.95CI.low = c(),
                            fisher.oddsratio.95CI.high = c(),
                            pearsons.phi.effectsize = c(),
                            pearsons.phi.effectsize.95CI.low = c(),
                            pearsons.phi.effectsize.95CI.high = c(),
                            sample.ID = c(),
                            bonferroni.corrected.chisq.test.pval = c(),
                            bonferroni.corrected.prop.test.pval = c(),
                            total.number.of.cancer.cells.in.sample = c(),
                            total.number.of.normal.cells.in.sample = c(),
                            total.number.of.cells.in.sample.after.low.quality.removed = c(),
                            stringsAsFactors = FALSE)
    old_table <- data.frame(Variant.probe.name = c(),#probe_table_list[[sample_ID]][ ,]),
                            Reference.probe.name = c(),
                            Number.of.variant.probes.in.cancer.cells = c(),
                            Number.of.reference.probes.in.cancer.cells = c(),
                            Number.of.variant.probes.in.non.cancer.cells = c(),
                            Number.of.reference.probes.in.non.cancer.cells = c(),
                            total.number.of.both.variant.site.probes.in.cancer.cells = c(),
                            total.number.of.both.variant.site.probes.in.normal.cells = c(),
                            proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = c(), # Proportion.test.value.1 # true positive rate
                            proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = c(), # Proportion.test.value.2
                            Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = c(),
                            Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = c(),
                            Chi.square.p.value = c(),
                            Proportion.test.95CI.low = c(),
                            Proportion.test.95CI.high = c(),
                            Proportion.test.p.value = c(),
                            fisher.oddsratio = c(),
                            fisher.oddsratio.95CI.low = c(),
                            fisher.oddsratio.95CI.high = c(),
                            pearsons.phi.effectsize = c(),
                            pearsons.phi.effectsize.95CI.low = c(),
                            pearsons.phi.effectsize.95CI.high = c(),
                            sample.ID = c(),
                            bonferroni.corrected.chisq.test.pval = c(),
                            bonferroni.corrected.prop.test.pval = c(),
                            total.number.of.cancer.cells.in.sample = c(),
                            total.number.of.normal.cells.in.sample = c(),
                            total.number.of.cells.in.sample.after.low.quality.removed = c(),
                            stringsAsFactors = FALSE)
    
    count = 0
    for (i in 1:length(input_table$sample_ID)) {
        sample_ID <- input_table$sample_ID[i]
        if (variant %in% probe_table_list[[sample_ID]][ , "Variant.probe.name" ]) {
            print(sample_ID)
            # find the index of that variant in the sample specific table
            sample_variant_index <- match(variant, probe_table_list[[sample_ID]][ , "Variant.probe.name" ])
            if (variant == "PIK3CA-p-E545K-ALT-21-A") {
                variant <- "PIK3CA-p-E545K-ALT-A"
                old_table <- read.table(paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"),sep='\t',header=T)
                rename = 1
                count = count + 1
                if (count == 1) {
                    old_table <- read.table(paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"),sep='\t',header=T)
                }
            } else if (variant == "KRAS-p-G12V-ALT-21-T") {
                variant <- "KRAS-p-G12V-ALT-A"
                old_table <- read.table(paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"),sep='\t',header=T)
                rename = 1
                count = count + 1
                if (count == 1) {
                    old_table <- read.table(paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"),sep='\t',header=T)
                }
            } else if (variant == "KRAS-p-G12D-ALT-21-A") {
                variant <- "KRAS-p-G12D-ALT-T"
                old_table <- read.table(paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"),sep='\t',header=T)
                rename = 1
                count = count + 1
                if (count == 1) {
                    old_table <- read.table(paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"),sep='\t',header=T)
                }
            }
            # generate a new index for that sample in the new table
            new_table_index = length(new_table$sample.ID) + 1
            # add all of the values for that sample.
            new_table[new_table_index, "Variant.probe.name"] <- variant
            new_table[new_table_index, "Reference.probe.name"] <- probe_table_list[[sample_ID]][sample_variant_index , "Reference.probe.name" ]
            new_table[new_table_index, "Chi.square.p.value"] <- probe_table_list[[sample_ID]][sample_variant_index , "Chi.square.p.value" ]
            new_table[new_table_index, "Proportion.test.p.value"] <- probe_table_list[[sample_ID]][sample_variant_index , "Proportion.test.p.value" ]
            new_table[new_table_index, "Proportion.test.95CI.low"] <- probe_table_list[[sample_ID]][sample_variant_index , "Proportion.test.95CI.low" ]
            new_table[new_table_index, "Proportion.test.95CI.high"] <- probe_table_list[[sample_ID]][sample_variant_index , "Proportion.test.95CI.high" ]
            new_table[new_table_index, "proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells" ]
            new_table[new_table_index, "proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells" ]
            new_table[new_table_index, "Number.of.variant.probes.in.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "Number.of.variant.probes.in.cancer.cells" ]
            new_table[new_table_index, "Number.of.reference.probes.in.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "Number.of.reference.probes.in.cancer.cells" ]
            new_table[new_table_index, "Number.of.variant.probes.in.non.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "Number.of.variant.probes.in.non.cancer.cells" ]
            new_table[new_table_index, "Number.of.reference.probes.in.non.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "Number.of.reference.probes.in.non.cancer.cells" ]
            new_table[new_table_index, "total.number.of.both.variant.site.probes.in.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "total.number.of.both.variant.site.probes.in.cancer.cells" ]
            new_table[new_table_index, "total.number.of.both.variant.site.probes.in.normal.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "total.number.of.both.variant.site.probes.in.normal.cells" ]
            new_table[new_table_index, "Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells" ]
            new_table[new_table_index, "Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells"] <- probe_table_list[[sample_ID]][sample_variant_index , "Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells" ]
            new_table[new_table_index, "total.number.of.cancer.cells.in.sample"] <- probe_table_list[[sample_ID]][sample_variant_index , "total.number.of.cancer.cells.in.sample" ]
            new_table[new_table_index, "total.number.of.normal.cells.in.sample"] <- probe_table_list[[sample_ID]][sample_variant_index , "total.number.of.normal.cells.in.sample" ]
            new_table[new_table_index, "total.number.of.cells.in.sample.after.low.quality.removed"] <- probe_table_list[[sample_ID]][sample_variant_index , "total.number.of.cells.in.sample.after.low.quality.removed" ]
            new_table[new_table_index, "bonferroni.corrected.chisq.test.pval"] <- probe_table_list[[sample_ID]][sample_variant_index , "bonferroni.corrected.chisq.test.pval" ]
            new_table[new_table_index, "bonferroni.corrected.prop.test.pval"] <- probe_table_list[[sample_ID]][sample_variant_index , "bonferroni.corrected.prop.test.pval" ]
            new_table[new_table_index, "fisher.oddsratio"] <- probe_table_list[[sample_ID]][sample_variant_index , "fisher.oddsratio" ]
            new_table[new_table_index, "fisher.oddsratio.95CI.low"] <- probe_table_list[[sample_ID]][sample_variant_index , "fisher.oddsratio.95CI.low" ]
            new_table[new_table_index, "fisher.oddsratio.95CI.high"] <- probe_table_list[[sample_ID]][sample_variant_index , "fisher.oddsratio.95CI.high" ]
            new_table[new_table_index, "pearsons.phi.effectsize"] <- probe_table_list[[sample_ID]][sample_variant_index , "pearsons.phi.effectsize" ]
            new_table[new_table_index, "pearsons.phi.effectsize.95CI.low"] <- probe_table_list[[sample_ID]][sample_variant_index , "pearsons.phi.effectsize.95CI.low" ]
            new_table[new_table_index, "pearsons.phi.effectsize.95CI.high"] <- probe_table_list[[sample_ID]][sample_variant_index , "pearsons.phi.effectsize.95CI.high" ]
            new_table[new_table_index, "sample.ID"] <- sample_ID
            if (rename == 1) {
                if (variant == "PIK3CA-p-E545K-ALT-A") {
                    #old_table <- rbind(old_table,new_table)
                    #write.table(old_table, file = paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"), sep="\t",quote=F)
                    variant <- "PIK3CA-p-E545K-ALT-21-A"
                    rename = 0
                } else if (variant == "KRAS-p-G12V-ALT-A") {
                    #write.table(old_table, file = paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"), sep="\t",quote=F)
                    variant <- "KRAS-p-G12V-ALT-21-T"
                    rename = 0
                } else if (variant == "KRAS-p-G12D-ALT-T") {
                    #write.table(old_table, file = paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"), sep="\t",quote=F)
                    variant <- "KRAS-p-G12D-ALT-21-A"
                    rename = 0
                }
            }
        }
    }
    if (variant == "PIK3CA-p-E545K-ALT-21-A") {
        variant <- "PIK3CA-p-E545K-ALT-A"
    } else if (variant == "KRAS-p-G12V-ALT-21-T") {
        variant <- "KRAS-p-G12V-ALT-A"
    } else if (variant == "KRAS-p-G12D-ALT-21-A") {
        variant <- "KRAS-p-G12D-ALT-T"
    }
    old_table <- rbind(old_table, new_table)
    # print("old_after_rbind")
    # print(old_table[,c("Variant.probe.name","Reference.probe.name","sample.ID")])
    write.table(old_table, file = paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"), sep="\t",quote=F)
    #write.table(new_table,file = paste0(out_dir,"/Variant_specific_results/",variant,"_probe_specificity_results_by_sample.tsv"), sep="\t",quote=F)
    merged_probe_table <- rbind(merged_probe_table,new_table)
}
write.table(merged_probe_table, file = paste0(out_dir,"/Variant_specific_results/","All_variants","_probe_specificity_results_by_sample.tsv"), sep="\t",quote=F)
