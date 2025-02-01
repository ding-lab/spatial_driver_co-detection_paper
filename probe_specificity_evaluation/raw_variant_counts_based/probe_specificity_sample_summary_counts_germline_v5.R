library(tidyverse)
library(epitools)
library(effectsize)
set.seed(1234)
args = commandArgs(trailingOnly=TRUE)
print(args)
input_table_path <- args[1]
output_table_name <- args[2]
all_variants <- read.table(input_table_path,sep='\t',header=T)
sample_mutation_table <- read.table("/diskmnt/Projects/Users/austins2/tools/xenium_snvs/sample_mutation_table_based_on_bulk_7.tsv",sep='\t',header=T,row.names = "X")
colnames(sample_mutation_table) <- c("variant","sample","germline")
mutation_vector <- sort(unique(sample_mutation_table$variant))
number_of_variants <- length(mutation_vector)
statistics_table_results <- data.frame(Variant.probe.name = mutation_vector,
                                       Reference.probe.name = rep(NA, number_of_variants),
                                       Number.of.variant.probes.in.cancer.cells = rep(NA, number_of_variants),
                                       Number.of.reference.probes.in.cancer.cells = rep(NA, number_of_variants),
                                       Number.of.variant.probes.in.non.cancer.cells = rep(NA, number_of_variants),
                                       Number.of.reference.probes.in.non.cancer.cells = rep(NA, number_of_variants),
                                       total.number.of.either.variant.site.probes.in.cancer.cells = rep(NA, number_of_variants),
                                       total.number.of.either.variant.site.probes.in.normal.cells = rep(NA, number_of_variants),
                                       true.positive.rate.restricted..proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = rep(NA, number_of_variants), # Proportion.test.value.1 # true positive rate
                                       false.positive.rate.restricted..proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = rep(NA, number_of_variants), # Proportion.test.value.2
                                       false.negative.rate.restricted..Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells = rep(NA, number_of_variants),
                                       true.negative.rate.restricted..Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells = rep(NA, number_of_variants),
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
                                       total.number.of.cancer.cells.in.sample.IDs.with.variant = rep(NA, number_of_variants),
                                       total.number.of.cancer.cells.in.non.variant.sample.IDs.and.normal.cells.across.all.samples = rep(NA, number_of_variants),
                                       total.number.of.cells.in.all.samples.after.low.quality.removed = rep(NA, number_of_variants),
                                       row.names = c(1:length(mutation_vector)))
for (i in 1:number_of_variants) {
    mutation = mutation_vector[i]
    print(mutation)
    variant_subset <- all_variants[(all_variants$Variant.probe.name == mutation),]
    Reference.probe.name <- unique(variant_subset$Reference.probe.name)[1]
    statistics_table_results[i, "Variant.probe.name"] <- mutation
    statistics_table_results[i, "Reference.probe.name"] <- Reference.probe.name
    sample_with_mutation_csv <- sample_mutation_table$sample[sample_mutation_table$variant == mutation]
    germline <- sample_mutation_table$germline[sample_mutation_table$variant == mutation]
    samples_with_mutation_vec <- unlist(strsplit(sample_with_mutation_csv, ","))
    print(samples_with_mutation_vec)
    not_samples_with_mutation_vec <- unique(variant_subset$sample.ID[(!(variant_subset$sample.ID %in% samples_with_mutation_vec))])
    not_samples_with_mutation <- paste(not_samples_with_mutation_vec, collapse=",")
    statistics_table_results[i, "sample.IDs.with.variant"] <- sample_with_mutation_csv
    statistics_table_results[i, "non.variant.sample.IDs"] <- not_samples_with_mutation
    variant_subset_cancer <- variant_subset[(variant_subset$sample.ID %in% samples_with_mutation_vec),]
    variant_subset_non_cancer <- variant_subset[(!(variant_subset$sample.ID %in% samples_with_mutation_vec)),]
    data <- data.frame( ALT = c(NA,NA), # Cancer_cell, Normal_Cell
                        WT = c(NA,NA), # Cancer_cell, Normal_Cell
                        row.names = c("Cancer_cell", "Normal_cell"))
    data["Cancer_cell","ALT"] = sum(variant_subset_cancer$Number.of.variant.probes.in.cancer.cells)
    data["Cancer_cell","WT"] = sum(variant_subset_cancer$Number.of.reference.probes.in.cancer.cells)
    data["Normal_cell","ALT"] = sum(variant_subset_cancer$Number.of.variant.probes.in.non.cancer.cells) + sum(variant_subset_non_cancer$Number.of.variant.probes.in.cancer.cells) + sum(variant_subset_non_cancer$Number.of.variant.probes.in.non.cancer.cells)
    data["Normal_cell","WT"] = sum(variant_subset_cancer$Number.of.reference.probes.in.non.cancer.cells) + sum(variant_subset_non_cancer$Number.of.reference.probes.in.cancer.cells) + sum(variant_subset_non_cancer$Number.of.reference.probes.in.non.cancer.cells)
    statistics_table_results[i, "total.number.of.cancer.cells.in.sample.IDs.with.variant"] <- sum(variant_subset_cancer$total.number.of.cancer.cells.in.sample)
    statistics_table_results[i, "total.number.of.cancer.cells.in.non.variant.sample.IDs.and.normal.cells.across.all.samples"] <- sum(variant_subset_cancer$total.number.of.normal.cells.in.sample) + sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
    statistics_table_results[i, "total.number.of.cells.in.all.samples.after.low.quality.removed"] <- sum(variant_subset_cancer$total.number.of.cells.in.sample.after.low.quality.removed) + sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
    if (germline == TRUE) {
        germline_data <- data 
        germline_data["Cancer_cell","ALT"] = sum(variant_subset_cancer$Number.of.variant.probes.in.cancer.cells) + sum(variant_subset_cancer$Number.of.variant.probes.in.non.cancer.cells)
        germline_data["Cancer_cell","WT"] = sum(variant_subset_cancer$Number.of.reference.probes.in.cancer.cells) + sum(variant_subset_cancer$Number.of.reference.probes.in.non.cancer.cells)
        germline_data["Normal_cell","ALT"] = sum(variant_subset_non_cancer$Number.of.variant.probes.in.cancer.cells) + sum(variant_subset_non_cancer$Number.of.variant.probes.in.non.cancer.cells)
        germline_data["Normal_cell","WT"] = sum(variant_subset_non_cancer$Number.of.reference.probes.in.cancer.cells) + sum(variant_subset_non_cancer$Number.of.reference.probes.in.non.cancer.cells)
        data <- germline_data
        statistics_table_results[i, "total.number.of.cancer.cells.in.sample.IDs.with.variant"] <- sum(variant_subset_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
        statistics_table_results[i, "total.number.of.cancer.cells.in.non.variant.sample.IDs.and.normal.cells.across.all.samples"] <- sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
        statistics_table_results[i, "total.number.of.cells.in.all.samples.after.low.quality.removed"] <- sum(variant_subset_cancer$total.number.of.cells.in.sample.after.low.quality.removed) + sum(variant_subset_non_cancer$total.number.of.cells.in.sample.after.low.quality.removed)
    }
    statistics_table_results[i, "Number.of.variant.probes.in.cancer.cells"] <- data["Cancer_cell","ALT"]
    statistics_table_results[i, "Number.of.reference.probes.in.cancer.cells"] <- data["Cancer_cell","WT"]
    statistics_table_results[i, "Number.of.variant.probes.in.non.cancer.cells"] <- data["Normal_cell","ALT"]
    statistics_table_results[i, "Number.of.reference.probes.in.non.cancer.cells"] <- data["Normal_cell","WT"]
    statistics_table_results[i, "total.number.of.either.variant.site.probes.in.cancer.cells"] <- data["Cancer_cell","ALT"] + data["Cancer_cell","WT"]
    statistics_table_results[i, "total.number.of.either.variant.site.probes.in.normal.cells"] <- data["Normal_cell","ALT"] + data["Normal_cell","WT"]
    if ((sum(data["Normal_cell",]) !=0) & (sum(data["Cancer_cell",]) !=0)) {
        print(data)
        variant_chisq_test1 <- chisq.test(data, simulate.p.value=TRUE, correct=FALSE)
        print("simulate.p.value=TRUE, correct=FALSE")
        print(variant_chisq_test1)
        x <- data[,'ALT']
        n <- rowSums(data)
        variant_prop_test <- prop.test(x, n , alternative = "greater", correct = FALSE) 
        print("prop_test alternative=greater")
        print(variant_prop_test)
        oddsratio.out = oddsratio(as.matrix(data), method = "fisher")
        print("Odds ratio")
        print(oddsratio.out)
        phi.effectsize <- phi(as.matrix(data), digits = 3)
        print("Effect size (Pearson's phi)")
        print(phi.effectsize)
        statistics_table_results[i, "true.positive.rate.restricted..proportion.test.value.1..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells"] <- variant_prop_test$estimate[1][[1]] # true.positive.rate if the sample actually has the mutation
        statistics_table_results[i, "false.positive.rate.restricted..proportion.test.value.2..Proportion.of.variant.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells"] <- variant_prop_test$estimate[2][[1]] # false.positive.rate if the sample actually has the mutation
        statistics_table_results[i, "false.negative.rate.restricted..Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.cancer.cells"] <- data["Cancer_cell","WT"] / (data["Cancer_cell","ALT"] + data["Cancer_cell","WT"]) #false.negative.rate if the sample actually has the mutation
        statistics_table_results[i, "true.negative.rate.restricted..Proportion.of.reference.allele.probes.out.of.total.number.of.variant.site.probes.in.normal.cells"] <- data["Normal_cell","WT"] / (data["Normal_cell","ALT"] + data["Normal_cell","WT"]) #true.negative.rate if the sample actually has the mutation
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
    }
}
write.table(statistics_table_results ,paste0("./",output_table_name),sep='\t',quote=FALSE)
