## This script extracts a cleaned count matrix from the 10X output and largely based on:
## https://github.com/MarioniLab/DropletUtils and
## https://bioconductor.org/packages/3.11/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html

#This script compiles the data for 5 x LMCs (LMC2B and LMC5-8 or HMC2B and HMC5-8) and 3 x NMCs (NMC5-7 or RB5-7)

library(DropletUtils)
#----------------#all#--------------------#
#make a list of path directories
samples <- c("../../data/Batch2/CellROutput/HM_053_SIGAD4/SIGAD4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/HM_024_SIGAE4/SIGAE4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/HM_057_SIGAF4/SIGAF4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/HM_058_SIGAG4/SIGAG4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/HM_021_SIGAH4/SIGAH4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/RB_P40_SIGAA4/SIGAA4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/RB_P43_SIGAB4/SIGAB4/outs/raw_feature_bc_matrix",
             "../../data/Batch2/CellROutput/RB_P45_SIGAC4/SIGAC4/outs/raw_feature_bc_matrix")
sample.names <- c("HMC2B","HMC5","HMC6","HMC7","HMC8","RB5","RB6","RB7")

#make single cell object
sce <- DropletUtils::read10xCounts(samples, sample.names=sample.names, col.names = TRUE)
colData(sce)$Barcode <- colnames(sce) <- paste(colData(sce)$Sample,colData(sce)$Barcode, sep="_")

#run empty drops
set.seed(100)
out <- emptyDrops(counts(sce))
is.cell <- out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)
#112854

table(Limited=out$Limited, Significant=is.cell)
#       Significant
#Limited  FALSE   TRUE
#  FALSE 213587   6985
#  TRUE       0 105869

sce

#formatting is.cell vector to be either FALSE or TRUE
is.cell[is.na(is.cell)] <- FALSE

#keeping the non-empty droplets
sce <- sce[,is.cell]

#keeping only the genes that show some expression
keep <- rowSums(counts(sce)) > 1
sce <- sce[keep,]
# dim of sce 26388 77271

#Find annotations for Ensembl gene ID's using the package "org.Hs.eg.db"
library(org.Hs.eg.db)
my.ids <- rownames(rowData(sce))
anno <- select(org.Hs.eg.db, keys=my.ids, keytype="ENSEMBL", column="SYMBOL")
anno <- anno[match(my.ids, anno$ENSEMBL),]
rowData(sce)$Symbol <- anno$SYMBOL

#add Mitochondrial gene annotations
mitoGenes <- read.table("../../data/miscData/MitoGenes.txt")
rowData(sce)$Mito <- rowData(sce)$ID %in% mitoGenes$V1

#add TF gene annotations
TFGenes <- read.table("../../data/miscData/TFcheckpoint_WithENSID.txt",sep="\t")
colnames(TFGenes) <- TFGenes[1,]
rowData(sce)$TF <- rowData(sce)$ID %in% TFGenes$ensembl_gene_id

#save data
cDat <- counts(sce)
pDat <- colData(sce)
fDat <- rowData(sce)

DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../../data/Batch2/ExpressionList.rds")
saveRDS(sce,file="../../data/Batch2/sce.rds")

sessionInfo()
