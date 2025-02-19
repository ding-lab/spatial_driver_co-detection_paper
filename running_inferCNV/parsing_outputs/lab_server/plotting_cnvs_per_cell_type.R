library(tidyverse)
library(optparse)
library(Seurat)
library(ggrastr)
library(ggpubr)
set.seed(1234)
print("libraries loaded properly")
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="path to cnv matrix that is output from infercnv_preprocessing_v3.py (use the one without '_T' in the name)",
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
    make_option(c("-r","--seurat_rds"),
                type="character",
                default=NA,
                help="seruat objects with cell type annotations",
                metavar="character"),
    make_option(c('-g','--gene_list'),
                type = 'character',
                default = "KRAS,TP53",
                help="a comma separated list of genes for which specific feature plots will be generated for each sample; e.g. 'GENE1,GENE2,GENE3,GENE4'",
                metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample=opt$sample_id
outPath=opt$output
tumor=opt$tumor
cnv_path=opt$input
rds_path=opt$seurat_rds
gene_csv = opt$gene_list
cell_type_path = opt$gene_list
obj <- readRDS(rds_path)
# the objects have the cell barcode format of sample_barcode (e.g. HT065B1-S1H1A3N1_AAACCCAAGCTCGTGC-1)
# they all also have the original_RNA_barcode column: (e.g. AAACCCAAGCTCGTGC-1)
cnv_table <- read.csv(cnv_path, sep="\t", header = TRUE, row.names = "index")
genes <- colnames(cnv_table)
# fixed so don't need to do this anymore with version 3.
# for (i in 1:ncol(cnv_table)) {
#     cnv_table[, i] <- gsub("x", "", as.character(cnv_table[, i]))
#     cnv_table[, i] <- as.numeric(cnv_table[, i])
# }
print("counting instances of CNVs")
cnv_in_rows <- function(x) {
    #x is a row of the dataframe
    #subset to only values greater than or equal to 2.0 this is y
    #y = x[x == 0]
    y = x[x < 1]
    #subset to only values equal to 0 this is z
    #z = x[x >= 2]
    z = x[x > 1]
    #count the total between y and z
    tally = length(y)+length(z)
}
Total = apply(cnv_table, 1, cnv_in_rows)
cnv_table <- cbind(cnv_table,Total_CNA = Total)
cnv_table$sample_barcode <- rownames(cnv_table)

# cell types are already in objects
print(head(obj$cell_type_xenium_variant_v1))
# adding cell types to metadata:
cnv_table <- cnv_table[order(row.names(cnv_table)), ]
print("showing part of cnv_table for QC purposes")
print(cnv_table[1:5,1:5])

cnv_table_small <- data.frame(cnv_total = cnv_table$Total_CNA, row.names = rownames(cnv_table))

obj <- AddMetaData(object = obj, metadata = cnv_table_small)
print(head(obj$cnv_total))

coordinates <- Embeddings(obj, reduction = "umap.30PC")
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length

print("plotting world domination")
pdf(paste(outPath,"/",sample,"_gene_CNV_Plots.pdf",sep=""),useDingbats=FALSE,height=7,width=8)
Idents(obj) <- "seurat_clusters"
print(rasterize(FeaturePlot(obj, features = 'cnv_total', raster = FALSE, order = TRUE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300))
Idents(obj) <- "cell_type_xenium_variant_v1"
VlnPlot(obj, features = 'cnv_total') + labs(title=paste(sample, "cnv_total on cell labels",sep = " "))
print(rasterize(DimPlot(object=obj, reduction = "umap.30PC", group.by = "cell_type_xenium_variant_v1",label=TRUE, label.size=6, repel = TRUE, raster = FALSE) + labs(title=paste(sample, "cell_type_xenium_variant_v1",sep = " ")) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300))

gene_list <- unlist(strsplit(gene_csv, ","))
if (gene_csv != "") {
    for (gene in gene_list) {
        print(gene)
        if (gene %in% colnames(cnv_table)) {
            if (!(paste(gene,"_inferCNV",sep="") %in% colnames(obj@meta.data))) {
                obj@meta.data[, paste(gene,'_inferCNV',sep='')] <- NA
            }
            gene_table <- data.frame(gene = cnv_table[, gene], row.names = rownames(cnv_table))
            print(head(gene_table))
            obj <- AddMetaData(object = obj, metadata = gene_table, col.name = paste(gene,'_inferCNV',sep=''))
            # use this one when working with the original individual objects
            p1 <- rasterize(FeaturePlot(obj, features = paste(gene,'_inferCNV',sep=''), raster = FALSE, order = TRUE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers="Point", dpi=300)
            # use this one when working with the subset objects
            # p1 <- FeaturePlot(obj, reduction = "umap.30PC", features = paste(gene,'_cnvs',sep=''), raster = FALSE, order = TRUE) + RotatedAxis()
            print(p1)
        }
    }
}
dev.off()
write.table(obj@meta.data,paste0(outPath,"/",sample,"_meta.data_with_cnvs.tsv"),sep='\t',quote=F)
dir.create(paste0(outPath,"/",sample,"_cell_type_by_gene_cnv_tables"))
if (gene_csv != "") {
    for (gene in gene_list) {
        if (gene %in% colnames(cnv_table)) {
            cell_type_cnv = as.data.frame.matrix(table(obj@meta.data[, paste(gene,'_inferCNV',sep='')], obj@meta.data[,"cell_type_xenium_variant_v1"]))
            print(cell_type_cnv)
            #cell_type_cnv = as.numeric(cell_type_cnv)
            cell_type_cnv$gene_sums = rowSums(cell_type_cnv, na.rm=TRUE)
            cell_type_cnv["cell_type_sums",] = colSums(cell_type_cnv, na.rm = T)
            write.table(cell_type_cnv, paste0(outPath,"/",sample,"_cell_type_by_gene_cnv_tables/",sample,"_",gene,"cell_type_by_cnv_table.tsv"),sep="\t",quote = F)

        }
    }
}
