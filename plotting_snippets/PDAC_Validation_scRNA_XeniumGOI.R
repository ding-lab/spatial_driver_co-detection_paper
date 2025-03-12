set.seed(123)
library(Matrix)
library(Seurat)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(viridis)
library(dplyr)
library(stringr)
library(gridExtra)
library(grid)
library(tidyverse)
library(dittoSeq)
library(pheatmap)
library(scCustomize)

all<-readRDS('Transition_Acinar_Ductal_Tumor_v3.rds')

newmetadata<-read.table("/diskmnt/Projects/HTAN_analysis/PDAC/PDAC/metadata/pdac_metadata_fullobj_withADMannotation.txt", row.names=1,sep = '\t', header = TRUE)
all <- AddMetaData(
  object = all,
  metadata = newmetadata)

#Update object
all<-UpdateSeuratObject(all)

##Subset Tumor Only
Idents(all)<-all@meta.data$cell_type_specific3
tumoronly<-subset(all,idents=c("PDAC"))

#Set idents to sample and subset only samples o finterest that have annotated
#G12V/G12D variants annotated
Idents(tumoronly)<-tumoronly@meta.data$sample_id
tumoronly_SOI<-subset(tumoronly,idents=c(
#"p.G12V",
"HT056P1_S1PA",
"HT056P1_S1PB",
"HT056P1_S1R1",
"HT060P1_S1PB",
"HT060P1_S1R1",
"HT125P1_S1H4",
"HT125P1_S1H8",
"HT185P1_S1H2",
"HT185P1_S1H3",
"HT185P1_S1H6",

#"p.G12D",   
"HT061P1_S1PA",
"HT061P1_S1PB",
"HT061P1_S1PC",
"HT061P1_S1R1",
"HT071P1_S1H1",
"HT071P1_S1H4",
"HT071P1_S1H5",
"HT071P1_S1H9",
"HT085P1_S1H2",
"HT085P1_S1H3",
"HT085P1_S1H4",
"HT122P1_S1H4",
"HT122P1_S1H5",
"HT122P1_S1H9",
"HT124P1_S1H1",
"HT124P1_S1H2",
"HT166P1_S1H2",
"HT166P1_S1H5",
"HT166P1_S1H6",
"HT191P1_S1H1",
"HT191P1_S1H2",
"HT191P1_S1H4",
"HT200P1_S1H2",
"HT200P1_S1H4"
))

tumoronly_SOI<-RenameIdents(tumoronly_SOI,
#"p.G12V",
"HT056P1_S1PA"="p.G12V",
"HT056P1_S1PB"="p.G12V",
"HT056P1_S1R1"="p.G12V",
"HT060P1_S1PB"="p.G12V",
"HT060P1_S1R1"="p.G12V",
"HT125P1_S1H4"="p.G12V",
"HT125P1_S1H8"="p.G12V",
"HT185P1_S1H2"="p.G12V",
"HT185P1_S1H3"="p.G12V",
"HT185P1_S1H6"="p.G12V",

#"p.G12D",   
"HT061P1_S1PA"="p.G12D",
"HT061P1_S1PB"="p.G12D",
"HT061P1_S1PC"="p.G12D",
"HT061P1_S1R1"="p.G12D",
"HT071P1_S1H1"="p.G12D",
"HT071P1_S1H4"="p.G12D",
"HT071P1_S1H5"="p.G12D",
"HT071P1_S1H9"="p.G12D",
"HT085P1_S1H2"="p.G12D",
"HT085P1_S1H3"="p.G12D",
"HT085P1_S1H4"="p.G12D",
"HT122P1_S1H4"="p.G12D",
"HT122P1_S1H5"="p.G12D",
"HT122P1_S1H9"="p.G12D",
"HT124P1_S1H1"="p.G12D",
"HT124P1_S1H2"="p.G12D",
"HT166P1_S1H2"="p.G12D",
"HT166P1_S1H5"="p.G12D",
"HT166P1_S1H6"="p.G12D",
"HT191P1_S1H1"="p.G12D",
"HT191P1_S1H2"="p.G12D",
"HT191P1_S1H4"="p.G12D",
"HT200P1_S1H2"="p.G12D",
"HT200P1_S1H4"="p.G12D")

tumoronly_SOI@meta.data$KRAS_Annotation<-Idents(tumoronly_SOI)

##update levels for plotting

tumoronly_SOI@meta.data$sample_id <- factor(tumoronly_SOI@meta.data$sample_id, 
                            levels=c(
#"p.G12V",
"HT056P1_S1PA",
"HT056P1_S1PB",
"HT056P1_S1R1",
"HT060P1_S1PB",
"HT060P1_S1R1",
"HT125P1_S1H4",
"HT125P1_S1H8",
"HT185P1_S1H2",
"HT185P1_S1H3",
"HT185P1_S1H6",

#"p.G12D",   
"HT061P1_S1PA",
"HT061P1_S1PB",
"HT061P1_S1PC",
"HT061P1_S1R1",
"HT071P1_S1H1",
"HT071P1_S1H4",
"HT071P1_S1H5",
"HT071P1_S1H9",
"HT085P1_S1H2",
"HT085P1_S1H3",
"HT085P1_S1H4",
"HT122P1_S1H4",
"HT122P1_S1H5",
"HT122P1_S1H9",
"HT124P1_S1H1",
"HT124P1_S1H2",
"HT166P1_S1H2",
"HT166P1_S1H5",
"HT166P1_S1H6",
"HT191P1_S1H1",
"HT191P1_S1H2",
"HT191P1_S1H4",
"HT200P1_S1H2",
"HT200P1_S1H4"))

tumoronly_SOI@meta.data$case <- factor(tumoronly_SOI@meta.data$case, 
                            levels=c(
#"p.G12V",
"HT056P1",
"HT060P1",
"HT125P1",
"HT185P1",

#"p.G12D",   
"HT061P1",
"HT071P1",
"HT085P1",
"HT122P1",
"HT124P1",
"HT166P1",
"HT191P1",
"HT200P1"))

##Gene List of Interest
xenium_GOI<-c("UBE2C","TOP2A","MKI67","CENPF","CDK1","CCNB2","VWA5A","MDM2","FAS","PCNA",
                                      "GPX2","SERPINB2","EDN1","CXCL6","CXCL2","LIF","IL7R","CD14",
                                      "CCR7","GPR183","RGS16","LAMP3","SELL","CAVIN1","MET","MYC","TFPI","CAV1",
                                      "VCAN","COL5A2","ANPEP","PMP22","TNC")

pdf("validatedGOI_scRNA.pdf",width=4, height=4,useDingbats=FALSE) # this is part of Figure 6D
DotPlot(tumoronly_SOI,assay="RNA",features=rev(validated_GOI),cols=("RdBu"),cluster.idents=TRUE)&coord_flip()
dev.off()

pdf("validatedGOI_scRNA_supp.pdf",width=9, height=5,useDingbats=FALSE) # this is part of Extended data figure 9b
DotPlot(tumoronly_SOI,assay="RNA",features=rev(validated_GOI),cols=("RdBu"),group.by="sample_id",cluster.idents=FALSE)&
    coord_flip()&theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

