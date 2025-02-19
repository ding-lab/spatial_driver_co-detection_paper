#!/usr/bin/env Rscript
library(Seurat)
library(optparse)
require(Matrix)
library(tidyverse)
library(dplyr)
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
obj <- subset(obj, subset = neoplasm_normal_unknown %in% c("neoplasm","normal")) #for this project the snRNAseq cell types and neoplasm vs normal labels are saved as part of the Seurat object rds file in the metadata
DefaultAssay(obj) <- 'RNA'
mat=GetAssayData(object = obj, assay = "RNA", slot = "counts")
mat=as.data.frame(mat)
write.table(mat,paste(outPath, "/raw_counts_matrix/", sample,'.raw_counts.tsv',sep=''),sep='\t',quote=FALSE)

obj$barcode <- rownames(obj@meta.data)
tmp_df <- obj@meta.data[,c("barcode","neoplasm_normal_unknown")]
write(paste(sample, "\t", "normal", "\t", tumor, sep=""), file=paste(outPath,"/annotations_file/",'reference_cells_3.txt',sep=""), append=TRUE)
write.table(tmp_df, paste(outPath,"/annotations_file/",sample,"/",sample,'.Barcode_Annotation.txt',sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
