library(tidyverse)
library(ggrastr)
library(ggpubr)
set.seed(1234)
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"
all_sample_summary_total_cancer <- rowSums(all_sample_summary[,c("Number.of.cancer.cells.with.the.variant.probe","Number.of.cancer.cells.with.the.reference.probe")], na.rm=TRUE)
all_sample_summary_total_normal <- rowSums(all_sample_summary[,c("Number.of.non.cancer.cells.with.the.variant.probe","Number.of.non.cancer.cells.with.the.reference.probe")], na.rm=TRUE)
all_sample_summary$total.cancer.cell.counts = all_sample_summary_total_cancer
all_sample_summary$total.normal.cell.counts = all_sample_summary_total_normal
all_sample_summary$prop.alt.in.cancer = all_sample_summary$Number.of.cancer.cells.with.the.variant.probe / all_sample_summary$total.cancer.cell.counts
all_sample_summary$prop.wt.in.cancer = all_sample_summary$Number.of.cancer.cells.with.the.reference.probe / all_sample_summary$total.cancer.cell.counts
all_sample_summary$prop.alt.in.non.cacner = all_sample_summary$Number.of.non.cancer.cells.with.the.variant.probe / all_sample_summary$total.normal.cell.counts
all_sample_summary$prop.wt.in.non.cacner = all_sample_summary$Number.of.non.cancer.cells.with.the.reference.probe / all_sample_summary$total.normal.cell.counts
iterables <- rownames(input_table)
all_cancer_alt_summary <- all_sample_summary[ , c("Variant.probe.name","prop.alt.in.cancer","sample.ID")]
colnames(all_cancer_alt_summary) <- c("probe_site","probe_proportion","sample.ID")
all_cancer_alt_summary$cell_type <- "Neoplastic cells"
all_cancer_alt_summary$probe_target <- "Alternate allele"
all_cancer_wt_summary <- all_sample_summary[ , c("Variant.probe.name","prop.wt.in.cancer","sample.ID")]
colnames(all_cancer_wt_summary) <- c("probe_site","probe_proportion","sample.ID")
all_cancer_wt_summary$cell_type <- "Neoplastic cells"
all_cancer_wt_summary$probe_target <- "Wildtype allele"
all_normal_alt_summary <- all_sample_summary[ , c("Variant.probe.name","prop.alt.in.non.cacner","sample.ID")]
colnames(all_normal_alt_summary) <- c("probe_site","probe_proportion","sample.ID")
all_normal_alt_summary$cell_type <- "Normal cells"
all_normal_alt_summary$probe_target <- "Alternate allele"
all_normal_wt_summary <- all_sample_summary[ , c("Variant.probe.name","prop.wt.in.non.cacner","sample.ID")]
colnames(all_normal_wt_summary) <- c("probe_site","probe_proportion","sample.ID")
all_normal_wt_summary$cell_type <- "Normal cells"
all_normal_wt_summary$probe_target <- "Wildtype allele"
# prop_complete_table <- rbind(all_cancer_alt_summary, all_normal_alt_summary)
prop_complete_table <- rbind(all_cancer_alt_summary, all_cancer_wt_summary, all_normal_alt_summary, all_normal_wt_summary)
prop_complete_table$cell_type <- factor(prop_complete_table$cell_type, levels = c("Normal cells","Neoplastic cells"))
prop_complete_table$probe_target <- factor(prop_complete_table$probe_target, levels = c("Wildtype allele","Alternate allele"))
write.table(paste0(out_dir,"/","prop_complete_table_v6_subclone.tsv"),sep='\t',quote=F)
library(ggrepel) # this contains the function geom_text_repel()
iterables = rownames(input_table)
for (i in iterables) {
    sample = input_table[i,"Sample_ID"]
    dir.create(paste0(out_dir,"/",sample))
    feature_csv = input_table[i,"feature_csv"]
    feature_list <- unlist(strsplit(feature_csv, ","))
    list_of_plots <- list()
    print(sample)
    print(feature_list)
    pdf(paste0(out_dir,"/",sample,"/",sample,"_box_plot_proportion_alt_cells_v6_subclone.pdf"),useDingbats = F, height=3, width = 6)
    for (k in 1:length(feature_list)) {
        probe <- feature_list[k] # the feature list only contains features that were detected in the sample by bulk WES mutation calling and manually validated by bam-rc
        print(probe)
        if (grepl("ALT", probe)) {
            print(probe)
            probe_list <- c(probe)
            if (probe == "KRAS-p-G12D-ALT-T") {
                probe_list <- c(probe_list,"KRAS-p-G12D-ALT-21-A")
            } else if (probe == "KRAS-p-G12D-ALT-21-A") {
                probe_list <- c(probe_list,"KRAS-p-G12D-ALT-T")
            } else if (probe == "KRAS-p-G12V-ALT-A") {
                probe_list <- c(probe_list,"KRAS-p-G12V-ALT-21-T")
            } else if (probe == "KRAS-p-G12V-ALT-21-T") {
                probe_list <- c(probe_list,"KRAS-p-G12V-ALT-A")
            }
            print(probe_list)
            mutant_samples <- c()
            for (p in probe_list) {
                print(p)
                mutant_samples <- c(mutant_samples,input_table$Sample_ID[grepl(p,input_table$mutant_csv)])
            }
            print(mutant_samples)
            probe_prop_table <- prop_complete_table[(prop_complete_table$probe_site %in% probe_list) & (prop_complete_table$probe_target == "Alternate allele"),] 
            probe_prop_table$mutant_sample <- rep("sample without mutation",dim(probe_prop_table)[1])
            probe_prop_table$mutant_sample[probe_prop_table$sample.ID %in% mutant_samples] <- "mutant sample"
            probe_prop_table <- probe_prop_table[(!(is.nan(probe_prop_table$probe_proportion))),]
            probe_prop_table$mutant_sample <- factor(probe_prop_table$mutant_sample, levels = c("mutant sample","sample without mutation"))
            p2 <- ggplot(probe_prop_table, aes(x = cell_type, y = probe_proportion)) +  
                geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = mutant_sample), position = position_jitter(seed = 1234)) +
                geom_text_repel(size = 3, min.segment.length = 0, show.legend = FALSE, force = 4, direction="y", max.overlaps = 100000, aes(segment.alpha = 1, x=cell_type, label=ifelse(sample.ID %in% mutant_samples, as.character(sample.ID),'')), position = position_jitter(seed = 1234)) +
                scale_y_continuous(limits = c(0,1.1), breaks = seq(from = 0, to = 1, by = 0.25)) +
                #geom_text(aes(x=probe_proportion, label=ifelse(sample.ID == sample, as.character(sample.ID),'')), vjust = -0.8, hjust = 0) +
                #geom_signif(comparisons=list(c("Normal cells", "Neoplastic cells")), annotations=paste0("adjusted p = ", formatC(all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"], format = "E", digits = 2)),
                #            y_position = 1.05, tip_length = 0, vjust=0) +
                theme_classic() +
                theme(text = element_text(color = "black", size = 8),
                      axis.text.x = element_text(color = "black", size = 8),
                      axis.text.y = element_text(color = "black", size = 8),
                      #legend.position="none",
                      legend.text = element_text(color = "black", size = 8),
                      axis.title=element_text(size=8))+
                theme(axis.ticks = element_line(color = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
                labs( y = paste0("Proportion"),
                      #x = "Treatment", 
                      title = paste0("Proportion of cells with variant probe"),
                      colour = probe) +#,
                #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Alternate.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
                theme(plot.title=element_text(size=8, hjust=0.25, face="bold", colour="black", vjust=1)) +
                theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black"))
            print(p2)
            list_of_plots[[k]] <- p2
        }
    }
    dev.off()
    print(length(list_of_plots))
}
