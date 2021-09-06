# Script to prepare cellranger data for downstream analysis

library(cellrangerRkit)
library(dplyr)
library(org.Hs.eg.db)

# ---- ReadData ----
# Read in CountMatrix_x.rds for NMC(RB) and LMC(HMC) and convert each to a normal matrix
gene_RB1_matrix <- readRDS("../../data/Batch3/RB_M39/1/CountMatrix_RB1.rds")
gene_RB2_matrix <- readRDS("../../data/Batch3/RB_M39/2/CountMatrix_RB2.rds")
gene_RB3_matrix <- readRDS("../../data/Batch3/RB_M39/3/CountMatrix_RB3.rds")
gene_HMC_matrix <- readRDS("../../data/Batch3/HMC_17025/CountMatrix_HMC.rds")
cDatRB1 <- as.matrix(gene_RB1_matrix)
cDatRB2 <- as.matrix(gene_RB2_matrix)
cDatRB3 <- as.matrix(gene_RB3_matrix)
cDatHMC <- as.matrix(gene_HMC_matrix)

#Checking row names (ENSEMBL id identical) in all matricies are equal
identical(rownames(cDatRB1),rownames(cDatRB2))
identical(rownames(cDatRB2),rownames(cDatRB3))
identical(rownames(cDatRB3),rownames(cDatHMC))

#Merge matricies into 1
cDat <- cbind(cDatRB1,cDatRB2,cDatRB3,cDatHMC)
#dim(cBind)=58735,25530

#reduce size of matrix (delete genes that have no expression by any cell)
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]
#after cleaning, dim(cDat)=32796,25530

# ---- Formatting ----

######fDat

#Find annotations for Ensembl gene ID's using the package "org.Hs.eg.db"
my.ids <- rownames(cDat)
anno <- select(org.Hs.eg.db, keys=my.ids, keytype="ENSEMBL", column="SYMBOL")
anno <- anno[match(my.ids, anno$ENSEMBL),]
head(anno)

#Generate the fDat data frame with the Ensembl gene id's and gene symbols
fDat <- data.frame(anno, stringsAsFactors = FALSE)
colnames(fDat) <- c('id', 'symbol')
rownames(fDat) <- fDat$id

#add Mitochondial gene annotations
mitoGenes <- read.table("../../data/MitoGenes.txt")
fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1

#######pDat

#I have annotated the data based on which cell type it is RB or HMC and from which chamber it was loaded onto the chip
pRB1 <- data.frame(colnames(cDatRB1), stringsAsFactors = FALSE)
colnames(pRB1) <- c('barcode')
pRB1$Sample <- c('RB1')
pRB1$State <- c('RB')

pRB2 <-data.frame(colnames(cDatRB2), stringsAsFactors = FALSE)
colnames(pRB2) <- c('barcode')
pRB2$Sample <- c('RB2')
pRB2$State <- c('RB')

pRB3 <-data.frame(colnames(cDatRB3), stringsAsFactors = FALSE)
colnames(pRB3) <- c('barcode')
pRB3$Sample <- c('RB3')
pRB3$State <- c('RB')

pHMC <-data.frame(colnames(cDatHMC), stringsAsFactors = FALSE)
colnames(pHMC) <- c('barcode')
pHMC$Sample <- c('HMC')
pHMC$State <- c('HMC')

#Adding all the dataframes together to make a super data frame...
pDat <- rbind(pRB1,pRB2,pRB3,pHMC,deparse.level=1)

####Save data
stopifnot(identical(rownames(fDat),rownames(cDat)) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../../data/Batch3/ExpressionList.rds")