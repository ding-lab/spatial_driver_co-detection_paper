library(tidyverse)
library(optparse)
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="path to tsv table with deg output of Seurat's FindAllMarkers function",
                metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_path = opt$input
if (grepl("/", input_path, fixed=TRUE)) {
    degs_file_name <- tail(unlist(strsplit(input_path,"[/]")), n=1)
    degs_file_name_no_ext <- tools::file_path_sans_ext(degs_file_name)
    output_dir <- dirname(input_path)
    
} else {
    output_dir <- "./"
    degs_file_name <- input_path
    degs_file_name_no_ext <- tools::file_path_sans_ext(degs_file_name)
}
print(degs_file_name)
print(degs_file_name_no_ext)
print(output_dir)
degs <- read.table(input_path, sep = '\t', header = T)
high_floor = 0.9
med_floor = 0.35
low_floor = 0
depleted_ceiling = 0
degs$shared_high_expression = ""
degs$shared_med_expression = ""
degs$shared_low_expression = ""
degs$depleted_expression = ""
genes.iterations <- rownames(degs)
for (gene.row in genes.iterations) {
    gene <- degs[gene.row,"gene"]
    current_cluster <- degs[gene.row,"cluster"]
    degs_subset <- degs[(degs$gene == gene),]
    num_cell_types <- length(degs_subset$cluster)
    if (num_cell_types > 1) {
        high_subset <- degs_subset[ ( degs_subset$avg_log2FC > high_floor ), "cluster" ]
        med_subset <- degs_subset[ ( high_floor > degs_subset$avg_log2FC ) & ( degs_subset$avg_log2FC > med_floor ), "cluster" ]
        low_subset <- degs_subset[ ( med_floor > degs_subset$avg_log2FC ) & ( degs_subset$avg_log2FC > low_floor ), "cluster" ]
        depleted_subset <- degs_subset[ ( depleted_ceiling > degs_subset$avg_log2FC ), "cluster" ]
        high_subset <- high_subset[ !high_subset %in% current_cluster ]
        med_subset <- med_subset[ !med_subset %in% current_cluster ]
        low_subset <- low_subset[ !low_subset %in% current_cluster ]
        depleted_subset <- depleted_subset[ !depleted_subset %in% current_cluster ]
        if( length(high_subset) != 0 ) {
            high_string <- paste(noquote(unique(high_subset)), collapse = ',')
            degs[gene.row,"shared_high_expression"] <- high_string
        }
        if( length(med_subset) != 0 ) {
            med_string <- paste(noquote(unique(med_subset)), collapse = ',')
            degs[gene.row,"shared_med_expression"] <- med_string
        }
        if( length(low_subset) != 0 ) {
            low_string <- paste(noquote(unique(low_subset)), collapse = ',')
            degs[gene.row,"shared_low_expression"] <- low_string
        }
        if( length(depleted_subset) != 0 ) {
            my_will_to_live_string <- paste(noquote(unique(depleted_subset)), collapse = ',')
            degs[gene.row,"depleted_expression"] <- my_will_to_live_string
        }
    }
}
write.table(degs, paste0(output_dir,'/',degs_file_name_no_ext,"_shared_markers_compared.tsv"), sep = '\t', quote = F)
