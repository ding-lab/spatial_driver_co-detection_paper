library(tidyverse)
library(epitools)
library(effectsize)
library(ggrepel)
library(ggrastr)
library(ggpubr)
set.seed(1234)
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table_path <- args[1]
output_file_name <- args[2] #"All_variants_by_sample_cohort_level_false_positive_rate.pdf"
sample_mutation_table_path <- args[3]
all_variants <- read.table(input_table_path,sep='\t',header=T)
sample_mutation_table <- read.table(sample_mutation_table_path,sep='\t',header=T,row.names = "X")
colnames(sample_mutation_table) <- c("variant","sample","germline")
mutation_vector <- sort(unique(sample_mutation_table$variant))
#plot order
mutation_vector <- c("ATM-p-Q628fs-ALT-T",'FANCA-p-E1240fs-ALT--','APC-p-S1298Ffs--3-ALT-T','BRCA2-p-A938fs-ALT--','ATM-p-Y1124---ALT-G','TP53-p-Y220C-ALT-C','CTNNB1-p-S45Y-ALT-A','DLGAP4-p-Q689---ALT-T','IDH1-p-R132C-ALT-A','IDH1-p-R132H-ALT-T','APC-p-R213---ALT-A','BAP1-p-T93A-ALT-C','DIS3-p-E429Q-ALT-G','DIS3-p-R780K-ALT-T','IRF4-p-Q300---ALT-T','NRAS-p-Q61H-ALT-G','TP53-p-R282W-ALT-A','BRAF-p-V600E-ALT-T','GNAS-p-R201C-ALT-T','KRAS-p-G12D-ALT-T','KRAS-p-G12V-ALT-A','KRAS-p-Q61H-ALT-G','PIK3CA-p-E545K-ALT-A','PELI1-p-R46T-ALT-G','SMG1-p-K3172T-ALT-G','EEF1A1-p-D442H-ALT-G','MAP1B-p-E611D-ALT-T','PGAP2-p---316Rext--40-ALT-C','TPD52L1-p-A21---ALT-T','LDHB-p-Q307---ALT-T')

number_of_variants <- length(mutation_vector)
dir.create("./plots/")
# statistics_table_results <- data.frame(Variant.probe.name = mutation_vector,
#                                        Reference.probe.name = rep(NA, number_of_variants),
#                                        Number.of.cancer.cells.with.the.variant.probe = rep(NA, number_of_variants),
#                                        Number.of.cancer.cells.with.the.reference.probe = rep(NA, number_of_variants),
#                                        Number.of.non.cancer.cells.with.the.variant.probe = rep(NA, number_of_variants),
#                                        Number.of.non.cancer.cells.with.the.reference.probe = rep(NA, number_of_variants),
#                                        total.number.of.cancer.cells.with.either.variant.site.probe = rep(NA, number_of_variants),
#                                        total.number.of.normal.cells.with.either.variant.site.probe = rep(NA, number_of_variants),
#                                        Chi.square.p.value = rep(NA, number_of_variants),
#                                        true.positive.rate.restricted..proportion.test.value.1..Proportion.of.cancer.cells.with.at.least.one.variant.allele.out.of.total.number.of.cancer.cells.with.variant.site.probe = rep(NA, number_of_variants),
#                                        false.positive.rate.restricted..proportion.test.value.2..Proportion.of.normal.cells.with.at.least.one.variant.allele.out.of.total.number.of.normal.cells.with.variant.site.probe = rep(NA, number_of_variants),
#                                        false.negative.rate.restricted..Proportion.of.cancer.cells.with.only.reference.allele.out.of.total.number.cancer.cells.with.variant.site.probe = rep(NA, number_of_variants),
#                                        true.negative.rate.restricted..Proportion.of.normal.cells.with.only.reference.allele.out.of.total.number.normal.cells.with.variant.site.probe = rep(NA, number_of_variants),
#                                        Proportion.of.cancer.cells.with.at.least.one.variant.allele.out.of.total.number.of.cancer.cells.with.and.without.variant.site.coverage..true.positive.rate.open = rep(NA, number_of_variants),
#                                        Proportion.of.normal.cells.with.at.least.one.variant.allele.out.of.total.number.of.normal.cells.with.and.without.variant.site.coverage..false.positive.rate.open = rep(NA, number_of_variants),
#                                        Proportion.test.95CI.low = rep(NA, number_of_variants),
#                                        Proportion.test.95CI.high = rep(NA, number_of_variants),
#                                        Proportion.test.p.value = rep(NA, number_of_variants),
#                                        fisher.oddsratio = rep(NA, number_of_variants),
#                                        fisher.oddsratio.95CI.low = rep(NA, number_of_variants),
#                                        fisher.oddsratio.95CI.high = rep(NA, number_of_variants),
#                                        pearsons.phi.effectsize = rep(NA, number_of_variants),
#                                        pearsons.phi.effectsize.95CI.low = rep(NA, number_of_variants),
#                                        pearsons.phi.effectsize.95CI.high = rep(NA, number_of_variants),
#                                        bonferroni.corrected.chisq.test.pval = rep(NA, number_of_variants),
#                                        bonferroni.corrected.prop.test.pval = rep(NA, number_of_variants),
#                                        sample.IDs.with.variant.validated.by.sequencing = rep(NA, number_of_variants), # sample.ID in all_variants table
#                                        non.variant.sample.IDs.validated.by.sequencing = rep(NA, number_of_variants), # sample.ID in all_variants table
#                                        total.number.of.cancer.cells.in.sample.IDs.with.variant = rep(NA, number_of_variants),
#                                        total.number.of.cancer.cells.in.non.variant.sample.IDs.and.normal.cells.across.all.samples = rep(NA, number_of_variants),
#                                        total.number.of.cells.in.all.samples.after.low.quality.removed = rep(NA, number_of_variants),
#                                        row.names = c(1:length(mutation_vector)))

list_of_plots = list()
#list_of_plots = vector('list', 30)
plot_counter = 1
weighted_true_positive_rate_mean_vec = c()
weighted_true_positive_rate_mean_vec = list()
number_of_samples_vector = c()
number_of_samples_vector = list()
for (i in 1:number_of_variants) {
    mutation = mutation_vector[i]
    print(mutation)
    variant_subset <- all_variants[(all_variants$Variant.probe.name == mutation),]
    Reference.probe.name <- unique(variant_subset$Reference.probe.name)[1]
    # statistics_table_results[i, "Variant.probe.name"] <- mutation
    # statistics_table_results[i, "Reference.probe.name"] <- Reference.probe.name
    sample_with_mutation_csv <- sample_mutation_table$sample[sample_mutation_table$variant == mutation]
    germline <- sample_mutation_table$germline[sample_mutation_table$variant == mutation]
    samples_with_mutation_vec <- unlist(strsplit(sample_with_mutation_csv, ","))
    print(samples_with_mutation_vec)
    not_samples_with_mutation_vec <- unique(variant_subset$sample.ID[(!(variant_subset$sample.ID %in% samples_with_mutation_vec))])
    not_samples_with_mutation <- paste(not_samples_with_mutation_vec, collapse=",")
    # statistics_table_results[i, "sample.IDs.with.variant"] <- sample_with_mutation_csv
    # statistics_table_results[i, "non.variant.sample.IDs"] <- not_samples_with_mutation
    variant_subset_cancer <- variant_subset[(variant_subset$sample.ID %in% samples_with_mutation_vec),]
    variant_subset_non_cancer <- variant_subset[(!(variant_subset$sample.ID %in% samples_with_mutation_vec)),]
    data <- data.frame( ALT = c(NA,NA), # Cancer_cell, Normal_Cell
                        WT = c(NA,NA), # Cancer_cell, Normal_Cell
                        row.names = c("Cancer_cell", "Normal_cell"))
    data["Cancer_cell","ALT"] = sum(variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe)
    data["Cancer_cell","WT"] = sum(variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe)
    data["Normal_cell","ALT"] = sum(variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe) + sum(variant_subset_non_cancer$Number.of.cancer.cells.with.the.variant.probe) + sum(variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.variant.probe)
    data["Normal_cell","WT"] = sum(variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe) + sum(variant_subset_non_cancer$Number.of.cancer.cells.with.the.reference.probe) + sum(variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
    # statistics_table_results[i, "total.number.of.cancer.cells.in.sample.IDs.with.variant"] <- sum(variant_subset_cancer$total.number.of.cancer.cells.in.sample)
    # statistics_table_results[i, "total.number.of.cancer.cells.in.non.variant.sample.IDs.and.normal.cells.across.all.samples"] <- sum(variant_subset_cancer$total.number.of.normal.cells.in.sample) + sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
    # statistics_table_results[i, "total.number.of.cells.in.all.samples.after.low.quality.removed"] <- sum(variant_subset_cancer$total.number.of.cells.in.sample.after.low.quality.removed) + sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
    if (germline) {
        germline_data <- data 
        germline_data["Cancer_cell","ALT"] = sum(variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe) + sum(variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe)
        germline_data["Cancer_cell","WT"] = sum(variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe) + sum(variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
        germline_data["Normal_cell","ALT"] = sum(variant_subset_non_cancer$Number.of.cancer.cells.with.the.variant.probe) + sum(variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.variant.probe)
        germline_data["Normal_cell","WT"] = sum(variant_subset_non_cancer$Number.of.cancer.cells.with.the.reference.probe) + sum(variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
        data <- germline_data
        # statistics_table_results[i, "total.number.of.cancer.cells.in.sample.IDs.with.variant"] <- sum(variant_subset_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
        # statistics_table_results[i, "total.number.of.cancer.cells.in.non.variant.sample.IDs.and.normal.cells.across.all.samples"] <- sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
        # statistics_table_results[i, "total.number.of.cells.in.all.samples.after.low.quality.removed"] <- sum(variant_subset_cancer$total.number.of.cells.in.sample.after.low.quality.removed) + sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
    }
    # statistics_table_results[i, "Number.of.cancer.cells.with.the.variant.probe"] <- data["Cancer_cell","ALT"]
    # statistics_table_results[i, "Number.of.cancer.cells.with.the.reference.probe"] <- data["Cancer_cell","WT"]
    # statistics_table_results[i, "Number.of.non.cancer.cells.with.the.variant.probe"] <- data["Normal_cell","ALT"]
    # statistics_table_results[i, "Number.of.non.cancer.cells.with.the.reference.probe"] <- data["Normal_cell","WT"]
    # statistics_table_results[i, "total.number.of.cancer.cells.with.either.variant.site.probe"] <- data["Cancer_cell","ALT"] + data["Cancer_cell","WT"]
    # statistics_table_results[i, "total.number.of.normal.cells.with.either.variant.site.probe"] <- data["Normal_cell","ALT"] + data["Normal_cell","WT"]
    if ((sum(data["Normal_cell",]) !=0) & (sum(data["Cancer_cell",]) !=0)) {
        print(data)
        variant_subset_cancer$true_positive_rate = (variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe) / (variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe)
        variant_subset_cancer$false_positive_rate = (variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe) / (variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
        variant_subset_cancer$false_negative_rate = (variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe) / (variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe)
        variant_subset_cancer$true_negative_rate = (variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe) / (variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
        variant_subset_cancer$false_positive_sample_coverage = "sufficient coverage"
        variant_subset_cancer$false_positive_sample_coverage[(variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe) < 50] = "low coverage"
        if (germline) {
            variant_subset_cancer$true_positive_rate = (variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe) / (variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
            variant_subset_cancer$false_positive_rate = NA
            variant_subset_cancer$false_negative_rate = (variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe) / (variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
            variant_subset_cancer$true_negative_rate = NA
        }
        
        variant_subset_non_cancer$true_positive_rate = NA
        variant_subset_non_cancer$false_positive_rate = (variant_subset_non_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.variant.probe) / (variant_subset_non_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.cancer.cells.with.the.reference.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
        variant_subset_non_cancer$false_negative_rate = NA
        variant_subset_non_cancer$true_negative_rate = (variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.reference.probe + variant_subset_non_cancer$Number.of.cancer.cells.with.the.reference.probe) / (variant_subset_non_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.cancer.cells.with.the.reference.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.reference.probe)
        variant_subset_non_cancer$false_positive_sample_coverage = "sufficient coverage"
        variant_subset_non_cancer$false_positive_sample_coverage[(variant_subset_non_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.cancer.cells.with.the.reference.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.variant.probe + variant_subset_non_cancer$Number.of.non.cancer.cells.with.the.reference.probe) < 50] = "low coverage"
        
        variant_subset_rates <- rbind(variant_subset_cancer, variant_subset_non_cancer)
        samples_with_mutation_vec
        variant_subset_rates$mutant = "sample without mutation detected by sequencing"
        variant_subset_rates$mutant[variant_subset_rates$sample.ID %in% samples_with_mutation_vec] = "sample with mutation detected by sequencing"
        
        variant_subset_rates$false_positive_sample_coverage_and_necrosis <- variant_subset_rates$false_positive_sample_coverage
        necrosis_sample_vector <- c("HT260C1-Th1K1L1U1","HT268B1-Th1H3L1U1","HT284P1-S1H1A1U1","WUPE82256U1-Fp1U1","WUPE25723U1-Fp1U1","SP001P1-Fp1U1","SN106H1-Ma1Fp2-7U1","SP002C1-Fp1U2") # Adjusted Negative Control Probe Rate greater than or equal to 0.7% or # Adjusted Negative Control Codeword Rate greater than or equal to 0.07% --- Look at the analysis_summary.html file or the metrics_summary.csv file and check the value for each sample
        variant_subset_rates$false_positive_sample_coverage_and_necrosis[variant_subset_rates$sample.ID %in% necrosis_sample_vector] <- paste0(variant_subset_rates$false_positive_sample_coverage_and_necrosis[variant_subset_rates$sample.ID %in% necrosis_sample_vector]," and necrotic")
        
        weighted_true_positive_rate_mean = sum(variant_subset_rates$true_positive_rate * ((variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe)/ sum(variant_subset_cancer$Number.of.cancer.cells.with.the.variant.probe + variant_subset_cancer$Number.of.cancer.cells.with.the.reference.probe, na.rm=T)), na.rm=T)
        weighted_true_positive_rate_mean_vec <- c(weighted_true_positive_rate_mean_vec, weighted_true_positive_rate_mean)
        if (i == 1) {
            output_df = variant_subset_rates
            weighted_true_positive_rate_mean_vec[[mutation]] = weighted_true_positive_rate_mean
        } else {
            output_df = rbind(output_df, variant_subset_rates)
            weighted_true_positive_rate_mean_vec[[mutation]] = weighted_true_positive_rate_mean
        }
        number_of_samples_vector <- c(number_of_samples_vector, length(unique(variant_subset_rates$sample.ID)))
        
        p1 <- ggplot(variant_subset_rates, aes(x = Variant.probe.name, y = false_positive_rate)) +
            geom_boxplot(outlier.shape = NA) + 
            geom_jitter(aes(colour = mutant, shape = false_positive_sample_coverage_and_necrosis), size =2, position = position_jitter(seed = 1234)) +
            geom_hline(yintercept = weighted_true_positive_rate_mean, linetype=2, linewidth=0.5, color = "#F3766E", show.legend = T) +
            scale_linetype_manual(name = " ", values = c(weighted_true_positive_rate_mean = 2), labels = c("weighted mean true positive rate")) +
            scale_shape_manual(values=c("sufficient coverage" = 16, "sufficient coverage and necrotic" = 1, "low coverage" = 17, "low coverage and necrotic" = 2))+
            scale_y_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, by = 0.2)) +
            theme_classic() +
            theme(text = element_text(color = "black", size = 8),
                  axis.text.x = element_text(color = "black", size = 8),
                  axis.text.y = element_text(color = "black", size = 8),
                  #legend.position="none",
                  legend.text = element_text(color = "black", size = 8),
                  axis.title=element_text(size=8))+
            theme(legend.position="bottom",
                  legend.text = element_text(color = "black", size = 8)) +
            theme(axis.ticks = element_line(color = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
            labs( y = paste0("Proportion"),
                  title = paste0(unique(variant_subset_rates$Variant.probe.name)),
                  colour = mutation,
                  subtitle = paste0("n = ", length(unique(variant_subset_rates$sample.ID)))) + 
            theme(plot.title=element_text(size=8, hjust=0.25, face="bold", colour="black", vjust=1)) +
            theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black"))
        
        list_of_plots[[plot_counter]] <- p1
        print(plot_counter)
        plot_counter = plot_counter + 1
        #list_of_plots = append(list_of_plots, p1)
        # if (i == 1) {
        #     print("looooooooooooooooop 1")
        #     list_of_plots[[plot_counter]] = p1
        #     plot_counter = plot_counter + 1
        # } else {
        #     print("bep\nbep\nbep\nbep\nbep\n")
        #     list_of_plots[[plot_counter]] = p1
        #     plot_counter = plot_counter + 1
            # list_of_plots[[ length(list_of_plots) + 1 ]] <- local({ggplot(variant_subset_rates, aes(x = Variant.probe.name, y = false_positive_rate)) +
            #         geom_boxplot(outlier.shape = NA) + 
            #         geom_jitter(aes(colour = mutant, shape = false_positive_sample_coverage_and_necrosis), size =2, position = position_jitter(seed = 1234)) +
            #         geom_hline(yintercept = weighted_true_positive_rate_mean, linetype=2, linewidth=0.5, color = "#F3766E", show.legend = T) +
            #         scale_linetype_manual(name = " ", values = c(weighted_true_positive_rate_mean = 2), labels = c("weighted mean true positive rate")) +
            #         scale_shape_manual(values=c("sufficient coverage" = 16, "sufficient coverage and necrotic" = 1, "low coverage" = 17, "low coverage and necrotic" = 2))+
            #         scale_y_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, by = 0.2)) +
            #         theme_classic() +
            #         theme(text = element_text(color = "black", size = 8),
            #               axis.text.x = element_text(color = "black", size = 8),
            #               axis.text.y = element_text(color = "black", size = 8),
            #               #legend.position="none",
            #               legend.text = element_text(color = "black", size = 8),
            #               axis.title=element_text(size=8))+
            #         theme(legend.position="bottom",
            #               legend.text = element_text(color = "black", size = 8)) +
            #         theme(axis.ticks = element_line(color = "black"),
            #               panel.grid.major = element_blank(),
            #               panel.grid.minor = element_blank(),
            #               panel.background = element_blank(),
            #               panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
            #         labs( y = paste0("Proportion"),
            #               title = paste0(unique(variant_subset_rates$Variant.probe.name)),
            #               colour = mutation,
            #               subtitle = paste0("n = ", length(unique(variant_subset_rates$sample.ID)))) + 
            #         theme(plot.title=element_text(size=8, hjust=0.25, face="bold", colour="black", vjust=1)) +
            #         theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black"))})
        # }
        pdf(paste0("./plots/",mutation, "_by_sample_cohort_level_false_positive_rate.pdf"), useDingbats = F, height = 3, width = 3)
        print(p1)
        dev.off()
    }
}

# plotting_order <- c("ATM-p-Q628fs-ALT-T",'FANCA-p-E1240fs-ALT--','APC-p-S1298Ffs--3-ALT-T','BRCA2-p-A938fs-ALT--','ATM-p-Y1124---ALT-G','TP53-p-Y220C-ALT-C','CTNNB1-p-S45Y-ALT-A','DLGAP4-p-Q689---ALT-T','IDH1-p-R132C-ALT-A','IDH1-p-R132H-ALT-T','APC-p-R213---ALT-A','BAP1-p-T93A-ALT-C','DIS3-p-E429Q-ALT-G','DIS3-p-R780K-ALT-T','IRF4-p-Q300---ALT-T','NRAS-p-Q61H-ALT-G','TP53-p-R282W-ALT-A','BRAF-p-V600E-ALT-T','GNAS-p-R201C-ALT-T','KRAS-p-G12D-ALT-T','KRAS-p-G12V-ALT-A','KRAS-p-Q61H-ALT-G','PIK3CA-p-E545K-ALT-A','PELI1-p-R46T-ALT-G','SMG1-p-K3172T-ALT-G','EEF1A1-p-D442H-ALT-G','MAP1B-p-E611D-ALT-T','PGAP2-p---316Rext--40-ALT-C','TPD52L1-p-A21---ALT-T','LDHB-p-Q307---ALT-T')
# num_plots = length(plotting_order)
# for (i in 1:num_plots) {
#     mutation = plotting_order[i]
#     variant_subset <- output_df[(output_df$Variant.probe.name == mutation),]
#     if (i == 1) {
#         output_df_2 = variant_subset
#     } else {
#         output_df_2 = rbind(output_df_2, variant_subset)
#     }
# }
# output_df = output_df_2
output_df$Variant.probe.name = factor(output_df$Variant.probe.name, levels = mutation_vector)

write.table(output_df, paste0("./plots/",output_file_name,"_data.tsv"),sep='\t',quote=F)
number_of_plots = length(list_of_plots)
print(length(number_of_plots))
pdf(paste0("./plots/",output_file_name), useDingbats = F, height = 3, width = 1*number_of_plots)
print(ggarrange(plotlist = list_of_plots, nrow = 1, ncol = number_of_plots))
dev.off()

p2 <- ggplot(output_df, aes(x = Variant.probe.name, y = false_positive_rate)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(colour = mutant, shape = false_positive_sample_coverage_and_necrosis), size =2, position = position_jitter(seed = 1234)) +
    #geom_hline(yintercept = weighted_true_positive_rate_mean, linetype=2, linewidth=0.5, color = "#F3766E", show.legend = T) +
    scale_shape_manual(values=c("sufficient coverage" = 16, "sufficient coverage and necrotic" = 1, "low coverage" = 17, "low coverage and necrotic" = 2)) +
    scale_y_continuous(limits = c(-0.1,1.1), breaks = seq(from = 0, to = 1, by = 0.2))# +
x_seg_start_vec = seq(0.75,number_of_plots-1+0.75,1)
x_seg_stop_vec = seq(1.25,number_of_plots+0.25,1)
y_seg_start_vec = weighted_true_positive_rate_mean_vec
y_seg_stop_vec = weighted_true_positive_rate_mean_vec
# label_x = 1:30
# label_y = rep(1, number_of_plots)
for (i in 1:number_of_plots) {
    mutation = mutation_vector[i]
    x_start = x_seg_start_vec[i]
    x_stop = x_seg_stop_vec[i]
    # y_start = y_seg_start_vec[[i]]
    # y_stop = y_seg_stop_vec[[i]]
    y_start = y_seg_start_vec[[mutation]]
    y_stop = y_seg_stop_vec[[mutation]]
    
    label_n = paste0("n = ",number_of_samples_vector[i])
    p2 <- p2 + geom_segment(x = x_start, xend = x_stop,  y = y_start, yend = y_stop, linetype=1, linewidth=0.5, color = "#F3766E", show.legend = F) +
        geom_text(label = label_n, x = i, y=1.05, color = "black", size = 2)
}
    #geom_segment(x = seq(0.75,number_of_plots-1+0.75,1), xend=seq(1.25,number_of_plots+0.25,1),  y = weighted_true_positive_rate_mean_vec, yend = weighted_true_positive_rate_mean_vec, linetype=2, linewidth=0.5, color = "#F3766E", show.legend = T) +
p2 <- p2 + scale_linetype_manual(name = " ", values = c(weighted_true_positive_rate_mean = 1), labels = c("weighted mean true positive rate")) +
    #geom_text(label = number_of_samples_vector, x = seq(from=1, to=number_of_plots, by=1), y=rep(1, number_of_plots), color = "black") +
    theme_classic() +
    theme(text = element_text(color = "black", size = 8),
          axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 8),
          #legend.position="none",
          legend.text = element_text(color = "black", size = 8),
          axis.title=element_text(size=8))+
    theme(legend.position="bottom",
          legend.text = element_text(color = "black", size = 8)) +
    theme(axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    labs( y = paste0("Proportion"),
          #title = paste0(unique(variant_subset_rates$Variant.probe.name)),
          colour = mutation) + #,
          #subtitle = paste0("n = ", length(unique(variant_subset_rates$sample.ID)))) + 
    theme(plot.title=element_text(size=8, hjust=0.25, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black"))
p2 <- rasterize(p2, layers="Point", dpi=300)
p2 <- rasterize(p2, layers="Polygon", dpi=300)
pdf(paste0("./plots/",output_file_name,"_one_plot.pdf"), useDingbats = F, height = 5, width = 12)
print(p2)
dev.off()

p3 <- ggplot(output_df, aes(x = Variant.probe.name, y = false_positive_rate)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(colour = mutant, shape = false_positive_sample_coverage_and_necrosis), size =2, position = position_jitter(seed = 1234)) +
    #geom_hline(yintercept = weighted_true_positive_rate_mean, linetype=2, linewidth=0.5, color = "#F3766E", show.legend = T) +
    scale_shape_manual(values=c("sufficient coverage" = 16, "sufficient coverage and necrotic" = 1, "low coverage" = 17, "low coverage and necrotic" = 2)) +
    scale_y_continuous(limits = c(-0.1,1.1), breaks = seq(from = 0, to = 1, by = 0.2))# +
x_seg_start_vec = seq(0.75,number_of_plots-1+0.75,1)
x_seg_stop_vec = seq(1.25,number_of_plots+0.25,1)
y_seg_start_vec = weighted_true_positive_rate_mean_vec
y_seg_stop_vec = weighted_true_positive_rate_mean_vec
# label_x = 1:30
# label_y = rep(1, number_of_plots)
# for (i in 1:number_of_plots) {
#     x_start = x_seg_start_vec[i]
#     x_stop = x_seg_stop_vec[i]
#     y_start = y_seg_start_vec[i]
#     y_stop = y_seg_stop_vec[i]
#     label_n = paste0("n = ",number_of_samples_vector[i])
#     p3 <- p3 + geom_segment(x = x_start, xend = x_stop,  y = y_start, yend = y_stop, linetype=1, linewidth=1, color = "#F3766E", show.legend = F) +
#         geom_text(label = label_n, x = i, y=1.05, color = "black", size = 2)
# }
for (i in 1:number_of_plots) {
    mutation = mutation_vector[i]
    x_start = x_seg_start_vec[i]
    x_stop = x_seg_stop_vec[i]
    # y_start = y_seg_start_vec[[i]]
    # y_stop = y_seg_stop_vec[[i]]
    y_start = y_seg_start_vec[[mutation]]
    y_stop = y_seg_stop_vec[[mutation]]
    
    label_n = paste0("n = ",number_of_samples_vector[i])
    p3 <- p3 + geom_segment(x = x_start, xend = x_stop,  y = y_start, yend = y_stop, linetype=1, linewidth=0.5, color = "#F3766E", show.legend = F) +
        geom_text(label = label_n, x = i, y=1.05, color = "black", size = 2)
}
for (i in 1:(number_of_plots-1)) {
    p3 <- p3 + geom_segment(x = i+0.5, xend = i+0.5,  y = -0.5, yend = 1.5, linetype=1, linewidth=0.5, color = "black", show.legend = F)
}
#geom_segment(x = seq(0.75,number_of_plots-1+0.75,1), xend=seq(1.25,number_of_plots+0.25,1),  y = weighted_true_positive_rate_mean_vec, yend = weighted_true_positive_rate_mean_vec, linetype=2, linewidth=0.5, color = "#F3766E", show.legend = T) +
p3 <- p3 + scale_linetype_manual(name = " ", values = c(weighted_true_positive_rate_mean = 1), labels = c("weighted mean true positive rate")) +
    #geom_text(label = number_of_samples_vector, x = seq(from=1, to=number_of_plots, by=1), y=rep(1, number_of_plots), color = "black") +
    theme_classic() +
    theme(text = element_text(color = "black", size = 8),
          axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 8),
          #legend.position="none",
          legend.text = element_text(color = "black", size = 8),
          axis.title=element_text(size=8))+
    theme(legend.position="bottom",
          legend.text = element_text(color = "black", size = 8)) +
    theme(axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    labs( y = paste0("Proportion"),
          #title = paste0(unique(variant_subset_rates$Variant.probe.name)),
          colour = mutation) + #,
    #subtitle = paste0("n = ", length(unique(variant_subset_rates$sample.ID)))) + 
    theme(plot.title=element_text(size=8, hjust=0.25, face="bold", colour="black", vjust=1)) +
    theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black"))

p3 <- rasterize(p3, layers="Point", dpi=300)
p3 <- rasterize(p3, layers="Polygon", dpi=300)
pdf(paste0("./plots/",output_file_name,"_one_plot_bars.pdf"), useDingbats = F, height = 5, width = 12)
print(p3)
dev.off()


