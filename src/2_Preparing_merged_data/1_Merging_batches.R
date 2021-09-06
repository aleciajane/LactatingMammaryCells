## ---- message=FALSE----------------------------------------------------------
#load libraries
library(plyr)
library(dplyr)
library(scater)
library(scran)
library(ggplot2)
library(knitr)

#specific for dendogram+heatmap
library(cowplot)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(umap)

library(batchelor)
library(scran)
library(scater)
library(org.Hs.eg.db)

#load data
sce3 <- readRDS("../data/Batch3/sce_final.rds")
colData(sce3)$Barcode <- colData(sce3)$barcode
colData(sce3)$Batch <- "B3"

sce1 <- readRDS("../data/Batch1/sce_final.rds")
colData(sce1)$Batch <- "B1"

sce2 <- readRDS("../data/Batch2/sce_final.rds")
colData(sce2)$Batch <- "B2"

inData <- read.csv("../data/lineage_markers.csv",stringsAsFactors = FALSE)
Genes <- inData$Gene

#add color scale
colour <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

## ----------------------------------------------------------------------------
sce3
sce1
sce2

## ----------------------------------------------------------------------------
universe <- Reduce(intersect, list(rownames(sce3),rownames(sce1), rownames(sce3)))
length(universe)

sce3 <-sce3[universe,]
sce1 <-sce1[universe,]
sce2 <-sce2[universe,]

## ----------------------------------------------------------------------------
rescaled <- multiBatchNorm(sce3,sce1,sce2)
sce3 <- rescaled[[1]]
sce1 <- rescaled[[2]]
sce2 <- rescaled[[3]]

## ----------------------------------------------------------------------------
sce2_HM <- subset(colData(sce2)$Barcode, colData(sce2)$State=="HM")
sce2_RB <- subset(colData(sce2)$Barcode, colData(sce2)$State=="RB")
sce2.HM <- sce2[,sce2_HM]
sce2.RB <- sce2[,sce2_RB]

sce1_HM <- subset(colData(sce1)$Barcode, colData(sce1)$State=="HM")
sce1_RB <- subset(colData(sce1)$Barcode, colData(sce1)$State=="RB")
sce1.HM <- sce1[,sce1_HM]
sce1.RB <- sce1[,sce1_RB]

sce3_HM <- subset(colData(sce3)$Barcode, colData(sce3)$State=="HM")
sce3_RB <- subset(colData(sce3)$Barcode, colData(sce3)$State=="RB")
sce3.HM <- sce3[,sce3_HM]
sce3.RB <- sce3[,sce3_RB]

## ----------------------------------------------------------------------------
dec1 <- modelGeneVar(sce3.RB)
dec2 <- modelGeneVar(sce1.RB)
dec3 <- modelGeneVar(sce2.RB)
dec4 <- modelGeneVar(sce3.HM)
dec5 <- modelGeneVar(sce1.HM)
dec6 <- modelGeneVar(sce2.HM)

combined.dec <- combineVar(dec1, dec2, dec3, dec4, dec5, dec6)
chosen.hvgs <-combined.dec$bio > 0
sum(chosen.hvgs)

## ----------------------------------------------------------------------------
rowData(sce3.RB) <- rowData(sce1.RB) <- rowData(sce2.RB) <- rowData(sce3.HM) <- rowData(sce1.HM) <- rowData(sce2.HM)

#formatting to include the same colnames for colData
keeps <- c("Sample","Barcode","State","sum","detected","subsets_Mito_sum","subsets_Mito_detected","subsets_Mito_percent","total","keep","sizeFactor","GraphClusters","Batch")   
colD1 <-colData(sce3.RB)[keeps]
colData(sce3.RB) <- colD1

colD2 <-colData(sce1.RB)[keeps]
colData(sce1.RB) <- colD2

colD3 <- colData(sce2.RB)[keeps]
colData(sce2.RB) <- colD3

colD4 <-colData(sce3.HM)[keeps]
colData(sce3.HM) <- colD4

colD5 <-colData(sce1.HM)[keeps]
colData(sce1.HM) <- colD5

colD5 <- colData(sce2.HM)[keeps]
colData(sce2.HM) <- colD5

#Removing reducedDim data from each experiment, so that they can be merged
reducedDim(sce3.RB, "PCA") <- NULL
reducedDim(sce3.RB, "UMAP") <- NULL
reducedDim(sce1.RB, "PCA") <- NULL
reducedDim(sce1.RB, "UMAP") <- NULL
reducedDim(sce2.RB, "PCA") <- NULL
reducedDim(sce2.RB, "UMAP") <- NULL
reducedDim(sce3.HM, "PCA") <- NULL
reducedDim(sce3.HM, "UMAP") <- NULL
reducedDim(sce1.HM, "PCA") <- NULL
reducedDim(sce1.HM, "UMAP") <- NULL
reducedDim(sce2.HM, "PCA") <- NULL
reducedDim(sce2.HM, "UMAP") <- NULL
sce3.RB <- removeAltExps(sce3.RB)
sce1.RB <- removeAltExps(sce1.RB)
sce2.RB <- removeAltExps(sce2.RB)
sce3.HM <- removeAltExps(sce3.HM)
sce1.HM <- removeAltExps(sce1.HM)
sce2.HM <- removeAltExps(sce2.HM)

## ----------------------------------------------------------------------------
set.seed(100)
all.sce <- cbind(sce2.RB,sce1.RB,sce3.RB,sce2.HM,sce1.HM,sce3.HM)

all.sce <- runPCA(all.sce,subset_row=chosen.hvgs, BSPARAM=BiocSingular::RandomParam())

#Determining rough number of clusters
snn.gr <- buildSNNGraph(all.sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters, Batch=all.sce$Batch)
tab

#Plotting the cells as a tSNE to see if batch effects are evident
all.sce <- runUMAP(all.sce, dimred="PCA", min_dist=0.5, n_neighbors=20)

#Add colours
sample.set <- c("HMC1", "HMC2", "HMC3","HMC4", "HMC2B","HMC5","HMC6","HMC7","HMC8", "HMC9","RB1", "RB2","RB3","RB4","RB5","RB6","RB7", "RB8")
colour.names <- c("salmon1", "peachpuff2", "burlywood", "#F8766D","peachpuff3", "salmon2", "pink1",  "indianred2", "indianred4", "hotpink2", "steelblue3", "slategray3", "paleturquoise3","#00BFC4","turquoise", "turquoise4", "lightblue1", "lightskyblue2")
batch.set <- c("B1", "B2", "B3")
colourbatch.names <- c("olivedrab4", "yellow3","navy")

plotReducedDim(all.sce, dimred="UMAP", colour_by="Sample") + scale_color_manual(values=setNames(colour.names, sample.set))
ggsave(filename="../data/PreCorrection_UMAP_Sample.pdf",width=8,height=7)
plotReducedDim(all.sce, dimred="UMAP", colour_by="Batch") + scale_color_manual(values=setNames(colourbatch.names, batch.set))
ggsave(filename="../data/PreCorrection_UMAP_Batch.pdf",width=8,height=7)

## ----------------------------------------------------------------------------
set.seed(100)

mnn.out <- fastMNN(sce2.RB,sce1.RB,sce3.RB,sce2.HM,sce1.HM,sce3.HM, d=50, k=20, subset.row=chosen.hvgs, correct.all=TRUE,
    BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

#preview of combined matrix
assay(mnn.out, "reconstructed")

rownames(mnn.out) <- rownames(all.sce)
colnames(mnn.out) <- colnames(all.sce)
counts(mnn.out, withDimnames=FALSE) <- counts(all.sce)
logcounts(mnn.out, withDimnames=FALSE) <- logcounts(all.sce)
rowData(mnn.out)$hvg <- chosen.hvgs

#Adding in gene symbol information for Ensembl gene ID's using the package "org.Hs.eg.db"
my.ids <- rownames(rowData(mnn.out))
anno <- select(org.Hs.eg.db, keys=my.ids, keytype="ENSEMBL", column="SYMBOL")
anno <- anno[match(my.ids, anno$ENSEMBL),]
rowData(mnn.out)$Symbol <- anno$SYMBOL

rowData(mnn.out)$has.symbol <- !is.na(rowData(mnn.out)$Symbol)
mnn.out <- mnn.out[rowData(mnn.out)$has.symbol,]
rownames(mnn.out) <- rowData(mnn.out)$Symbol

#identical(colnames(mnn.out),colnames(all.sce))
colData(mnn.out)$Barcode <- colData(all.sce)$Barcode
colData(mnn.out)$Sample <- colData(all.sce)$Sample
colData(mnn.out)$State <- colData(all.sce)$State
colData(mnn.out)$sum <- colData(all.sce)$sum
colData(mnn.out)$detected <- colData(all.sce)$detected
colData(mnn.out)$subsets_Mito_percent <- colData(all.sce)$subsets_Mito_percent
colData(mnn.out)$Batches <- colData(all.sce)$Batch
colData(mnn.out)$OldClust <- paste0(colData(all.sce)$Batch,"-",colData(all.sce)$GraphClusters)

## ----------------------------------------------------------------------------
#Clustering the samples and verifying that the resulting clusters are contributed by multiple batches 
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_louvain(snn.gr)$membership
colData(mnn.out)$Clusters <- factor(paste0("C",clusters.mnn))
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn

#Now examining how variable the clusters are, so we can diagnose whether there are many batch specific clusters
# Avoid minor difficulties with the 'table' class.
tab.mnn <- unclass(tab.mnn)
# Using a large pseudo.count to avoid unnecessarily large variances when the counts are low.
norm <- normalizeCounts(tab.mnn, pseudo.count=10)
# Ranking clusters by the largest variances.
rv <- rowVars(norm)
DataFrame(Batch=tab.mnn, var=rv)[order(rv, decreasing=TRUE),]

#-------#Now generating UMAP plots to examine the integration#-------#
set.seed(100)
mnn.out <- runPCA(mnn.out, dimred="corrected")
mnn.out <- runUMAP(mnn.out, dimred="PCA", min_dist=0.5, n_neighbors=20)
snn.gr <- buildSNNGraph(mnn.out, use.dimred="PCA")
clusters.mnn <- igraph::cluster_louvain(snn.gr)$membership
colData(mnn.out)$GraphClusters <- factor(paste0("C",clusters.mnn))
mnn.out$batch <- factor(mnn.out$batch)

PCA <- as.data.frame(reducedDim(mnn.out, "PCA"))
colData(mnn.out)$PCA1 <- PCA$PC1
colData(mnn.out)$PCA2 <- PCA$PC2
UMAP <- as.data.frame(reducedDim(mnn.out, "UMAP"))
colData(mnn.out)$UMAP1 <- UMAP$V1
colData(mnn.out)$UMAP2 <- UMAP$V2

#PCA
plotReducedDim(mnn.out, dimred="PCA", colour_by="Batches") + scale_color_manual(values=setNames(colourbatch.names, batch.set))
ggsave(filename="../data/PostCorrect_PCA_batch.pdf",width=8,height=7) 

#UMAPs coloured by batch, graphclusters and key marker genes
plotReducedDim(mnn.out, dimred="UMAP", colour_by="Batches") + scale_color_manual(values=setNames(colourbatch.names, batch.set))
ggsave(filename="../data/PostCorrect_UMAP_batch.pdf",width=8,height=7) 

plotReducedDim(mnn.out, dimred="UMAP", colour_by="Clusters")
ggsave(filename="../data/PostCorrect_UMAP_Clusters.pdf",width=8,height=7)

plotReducedDim(mnn.out, dimred="UMAP", colour_by="GraphClusters")
ggsave(filename="../data/PostCorrect_UMAP_Graphclusters.pdf",width=8,height=7)

plotReducedDim(mnn.out, dimred="UMAP", colour_by="Sample") + scale_color_manual(values=setNames(colour.names, sample.set))
ggsave(filename="../data/PostCorrect_UMAP_Sample.pdf",width=8,height=7)

#umap coloured by genes
myGraphs <- list()
pdf("../data/PostCorrect_UMAP_marker_genes.pdf", width=8, height=6.5)
for (Gene in Genes) {
    myGene <- as.numeric(logcounts(mnn.out)[Gene,])
    tsne.gene <- ggcells(mnn.out, aes(x=UMAP1, y=UMAP2, color=myGene))+
        scale_color_gradientn(colors=colour) +
        geom_point() +
        ggtitle(Gene) +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_blank(), legend.position = c(0.9, 0.9), legend.title = element_blank(), plot.title = element_text(family="sans", face="bold.italic", size=45, vjust=-5, hjust=0.1), axis.text.x= element_blank(), axis.text.y = element_blank(), axis.title.x= element_blank(), axis.title.y= element_blank(), axis.ticks = element_blank())
    plot(tsne.gene)
    myGraphs[[Gene]] <- tsne.gene
}
dev.off()

#Now checking the proportion of variance within each batch that is lost during MNN correction (i.e. within-batch variance that is removed during orthogonalization) If there is a large proportion of variance (>10%), this suggests that the correction is removing genuine biological variation.
metadata(mnn.out)$merge.info$lost.var

## ----------------------------------------------------------------------------
sessionInfo()

#Saving data
saveRDS(mnn.out,file="../data/sce_all_final_MilkCorrected_2.rds")
