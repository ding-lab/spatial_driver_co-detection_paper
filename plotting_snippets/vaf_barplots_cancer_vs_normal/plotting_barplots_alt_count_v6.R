library(tidyverse)
library(ggrastr)
library(ggpubr)
library(scales)
set.seed(1234)
print("loading inputs")
all_sample_summary <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/counts_based/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep='\t',header=T)
all_sample_summary_cell <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/probe_specificity_test/new/removed_unknown_v6_subclone/Variant_specific_results/All_variants_probe_specificity_results_by_sample.tsv",sep="\t",header=T)
input_table <- read.table("/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/neo_norm_unk_table_v7.tsv",sep='\t',header=T)
out_dir <- "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/cell_types/neoplastic_normal_unknown/coord_test/darkblue_ceiling1_v7/"

# adding in total cells where varaint site was detected to plot
print("adding count of cells where variant was detected")
all_samples = unique(all_sample_summary_cell$sample.ID)
all_features = unique(all_sample_summary_cell$Variant.probe.name)
all_sample_summary$total.cancer.cells.where.variant.site.was.sampled = NA
all_sample_summary$total.non.cancer.cells.where.variant.site.was.sampled = NA
for (sample_ID in all_samples) {
    for (feature in all_features) {
        total_cancer_with_variant_site_detected = all_sample_summary_cell$total.number.of.cancer.cells.with.either.variant.site.probe[(all_sample_summary_cell$sample.ID == sample_ID) & (all_sample_summary_cell$Variant.probe.name == feature)]
        all_sample_summary$total.cancer.cells.where.variant.site.was.sampled[(all_sample_summary$sample.ID == sample_ID) & (all_sample_summary$Variant.probe.name == feature)] = total_cancer_with_variant_site_detected
        
        total_normal_with_variant_site_detected = all_sample_summary_cell$total.number.of.normal.cells.with.either.variant.site.probe[(all_sample_summary_cell$sample.ID == sample_ID) & (all_sample_summary_cell$Variant.probe.name == feature)]
        all_sample_summary$total.non.cancer.cells.where.variant.site.was.sampled[(all_sample_summary$sample.ID == sample_ID) & (all_sample_summary$Variant.probe.name == feature)] = total_normal_with_variant_site_detected
    }
}

print("re-adding proportions to include those not previously included")
all_sample_summary$prop.var.in.cancer = all_sample_summary$Number.of.variant.probes.in.cancer.cells / all_sample_summary$total.number.of.both.variant.site.probes.in.cancer.cells
all_sample_summary$prop.ref.in.cancer = all_sample_summary$Number.of.reference.probes.in.cancer.cells / all_sample_summary$total.number.of.both.variant.site.probes.in.cancer.cells
all_sample_summary$prop.var.in.non.cacner = all_sample_summary$Number.of.variant.probes.in.non.cancer.cells / all_sample_summary$total.number.of.both.variant.site.probes.in.normal.cells
all_sample_summary$prop.ref.in.non.cacner = all_sample_summary$Number.of.reference.probes.in.non.cancer.cells / all_sample_summary$total.number.of.both.variant.site.probes.in.normal.cells
iterables <- rownames(input_table)

print("looping over samples for plotting each mutant variant")
for (i in iterables) {
    sample = input_table[i,"Sample_ID"]
    print(sample)
    dir.create(paste0(out_dir,"/",sample))
    feature_csv = input_table[i,"feature_csv_table_plots"]
    feature_list <- unlist(strsplit(feature_csv, ","))
    print(feature_list)
    # converting input counts to proportions for input to ggplot
    sample_summary <- all_sample_summary[c(all_sample_summary$sample.ID == sample),]
    cancer_alt_sample_summary <- sample_summary[c(sample_summary$sample.ID == sample),c("Variant.probe.name","prop.var.in.cancer","sample.ID","Number.of.variant.probes.in.cancer.cells","total.cancer.cells.where.variant.site.was.sampled")]
    colnames(cancer_alt_sample_summary) <- c("probe_site","probe_proportion","sample.ID","probe_count","cell_count")
    cancer_alt_sample_summary$cell_type <- "Cancer cells"
    cancer_alt_sample_summary$probe_target <- "Variant allele"
    cancer_wt_sample_summary <- sample_summary[c(sample_summary$sample.ID == sample),c("Variant.probe.name","prop.ref.in.cancer","sample.ID","Number.of.reference.probes.in.cancer.cells","total.cancer.cells.where.variant.site.was.sampled")]
    colnames(cancer_wt_sample_summary) <- c("probe_site","probe_proportion","sample.ID","probe_count","cell_count")
    cancer_wt_sample_summary$cell_type <- "Cancer cells"
    cancer_wt_sample_summary$probe_target <- "Reference allele"
    normal_alt_sample_summary <- sample_summary[c(sample_summary$sample.ID == sample),c("Variant.probe.name","prop.var.in.non.cacner","sample.ID","Number.of.variant.probes.in.non.cancer.cells","total.non.cancer.cells.where.variant.site.was.sampled")]
    colnames(normal_alt_sample_summary) <- c("probe_site","probe_proportion","sample.ID","probe_count","cell_count")
    normal_alt_sample_summary$cell_type <- "Normal cells"
    normal_alt_sample_summary$probe_target <- "Variant allele"
    normal_wt_sample_summary <- sample_summary[c(sample_summary$sample.ID == sample),c("Variant.probe.name","prop.ref.in.non.cacner","sample.ID","Number.of.reference.probes.in.non.cancer.cells","total.non.cancer.cells.where.variant.site.was.sampled")]
    colnames(normal_wt_sample_summary) <- c("probe_site","probe_proportion","sample.ID","probe_count","cell_count")
    normal_wt_sample_summary$cell_type <- "Normal cells"
    normal_wt_sample_summary$probe_target <- "Reference allele"
    sample_complete_table <- rbind(cancer_alt_sample_summary, cancer_wt_sample_summary, normal_alt_sample_summary, normal_wt_sample_summary)
    sample_complete_table$cell_count = paste0("n_cell=",sample_complete_table$cell_count)
    sample_complete_table$cell_type <- factor(sample_complete_table$cell_type, levels = c("Normal cells","Cancer cells"))
    sample_complete_table$probe_target <- factor(sample_complete_table$probe_target, levels = c("Variant allele","Reference allele"))
    # plotting bar plot of proportions (could do stacked but I think staggered looks better because then I can have the pvalue as part of the plot instead of as the subtitle
    list_of_plots <- list()
    for (k in 1:length(feature_list)) {
        alt_probe <- feature_list[k]
        print(alt_probe)
        if (grepl("ALT",alt_probe)) {
            sample_probe_table_complete <- sample_complete_table[c(sample_complete_table$probe_site == alt_probe),]
            total_probe_in_normal <- sum(sample_probe_table_complete$probe_count[(sample_probe_table_complete$cell_type == "Normal cells")])
            total_probe_in_cancer <- sum(sample_probe_table_complete$probe_count[(sample_probe_table_complete$cell_type == "Cancer cells")])
            sample_probe_table_complete$probe_total <- NA
            sample_probe_table_complete$probe_total[(sample_probe_table_complete$cell_type == "Cancer cells")] <- total_probe_in_cancer
            sample_probe_table_complete$probe_total[(sample_probe_table_complete$cell_type == "Normal cells")] <- total_probe_in_normal
            print(sample_probe_table_complete)
            # df_max_round = round(max(sample_probe_table_complete$probe_total), -2)+100
            df_max_round = max(sample_probe_table_complete$probe_total)
            #plot_max = 1.2*df_max_round
            cols <- c("Reference allele" = hue_pal()(length(unique(sample_probe_table_complete$probe_target)))[2], "Variant allele" = hue_pal()(length(unique(sample_probe_table_complete$probe_target)))[1])
            p1 <- ggplot(sample_probe_table_complete, aes(x = cell_type, y = probe_count, fill = probe_target)) + 
                geom_bar(stat = "identity", position = "stack") +
                scale_y_continuous(limits = c(0,1.2*df_max_round)) +
                scale_fill_manual(values = cols) +
                # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
                # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
                geom_signif(comparisons=list(c("Normal cells", "Cancer cells")), annotations=paste0("Adj. prop. test p = ", formatC(all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Variant.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"], format = "E", digits = 2)),
                            y_position = 1.1*df_max_round, tip_length = 0, vjust=0) +
                geom_text(data = subset(sample_probe_table_complete, probe_target == "Variant allele"),
                          aes(label = cell_count, x = cell_type, y = probe_total),
                          color = "black", nudge_y = max(sample_probe_table_complete$probe_total)*0.1) +
                theme_classic() +
                theme(axis.text.x = element_text(color = "black", size = 12),
                      axis.text.y = element_text(color = "black", size = 12),
                      #legend.position="none",
                      axis.title=element_text(size=12),
                      axis.ticks = element_line(color = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
                labs( y = paste0("Probe count"),
                      #x = "Treatment", 
                      title = paste0("Count of probes in cells"),
                      fill = alt_probe) +#,
                #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Variant.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
                theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
                theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
            list_of_plots[[k]] <- p1
        }
    }
    pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_count_based_counts_p_value_alt_v6.pdf"),useDingbats = F, height=5, width = 5)
    for ( k in 1:length(list_of_plots)) {
        print(list_of_plots[[k]])
    }
    dev.off()
    list_of_plots <- list()
    for (k in 1:length(feature_list)) {
        alt_probe <- feature_list[k]
        print(alt_probe)
        if (grepl("ALT",alt_probe)) {
            sample_probe_table_complete <- sample_complete_table[c(sample_complete_table$probe_site == alt_probe),]
            total_probe_in_normal <- sum(sample_probe_table_complete$probe_count[(sample_probe_table_complete$cell_type == "Normal cells")])
            total_probe_in_cancer <- sum(sample_probe_table_complete$probe_count[(sample_probe_table_complete$cell_type == "Cancer cells")])
            sample_probe_table_complete$probe_total <- NA
            sample_probe_table_complete$probe_total[(sample_probe_table_complete$cell_type == "Cancer cells")] <- total_probe_in_cancer
            sample_probe_table_complete$probe_total[(sample_probe_table_complete$cell_type == "Normal cells")] <- total_probe_in_normal
            print(sample_probe_table_complete)
            df_max_round = 1
            #plot_max = 1.2*df_max_round
            cols <- c("Reference allele" = hue_pal()(length(unique(sample_probe_table_complete$probe_target)))[2], "Variant allele" = hue_pal()(length(unique(sample_probe_table_complete$probe_target)))[1])
            p1 <- ggplot(sample_probe_table_complete, aes(x = cell_type, y = probe_proportion, fill = probe_target)) +
                geom_bar(stat = "identity", position = "stack") +
                scale_y_continuous(limits = c(0,1.2*df_max_round)) +
                scale_fill_manual(values = cols) +
                # how to format decimals in scientific notation: https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-labe-at-most-two-decimals
                # how to add significance bars to plot: https://stackoverflow.com/questions/17084566/put-stars-on-ggplot-barplots-and-boxplots-to-indicate-the-level-of-significanc
                geom_signif(comparisons=list(c("Normal cells", "Cancer cells")), annotations=paste0("Adj. prop. test p = ", formatC(all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Variant.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"], format = "E", digits = 2)),
                            y_position = 1.1*df_max_round, tip_length = 0, vjust=0) +
                geom_text(data = subset(sample_probe_table_complete, probe_target == "Variant allele"),
                          aes(label = cell_count, x = cell_type, y = df_max_round),
                          color = "black", nudge_y = df_max_round*0.1) +
                theme_classic() +
                theme(axis.text.x = element_text(color = "black", size = 12),
                      axis.text.y = element_text(color = "black", size = 12),
                      #legend.position="none",
                      axis.title=element_text(size=12),
                      axis.ticks = element_line(color = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) +
                labs( y = paste0("Proportion"),
                      #x = "Treatment",
                      title = paste0("Proportion of probes in cells"),
                      fill = alt_probe) +#,
                #subtitle= paste0("Proportion test adjusted p = ",all_sample_summary[(all_sample_summary$sample.ID == sample)& (all_sample_summary$Variant.probe.name == alt_probe),"bonferroni.corrected.prop.test.pval"])) +
                theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=1)) +
                theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black"))
            list_of_plots[[k]] <- p1
        }
    }
    pdf(paste0(out_dir,"/",sample,"/",sample,"_barplots_count_based_proportions_p_value_alt_v6.pdf"),useDingbats = F, height=5, width = 5)
    for ( k in 1:length(list_of_plots)) {
        print(list_of_plots[[k]])
    }
    dev.off()
}
