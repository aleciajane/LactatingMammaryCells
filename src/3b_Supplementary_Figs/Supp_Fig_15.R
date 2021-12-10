## ---- message=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(plyr)
library(dplyr)
library(diffusionMap)

#Bioc Packages
library(scran)
library(scater)
library(destiny)
library(RColorBrewer)

##load data
sce <- readRDS("../../data/sce_all_nospike_2.rds")
ClustCol <- as.data.frame(read.csv("../../data/Clusters_noSpike_2.csv",stringsAsFactors = FALSE))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sce


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#UMAP classes colours
classes <- ClustCol$Identity
colours.classes <- ClustCol$Identity.colours

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sce.sub <- sce[,(colData(sce)$Identity=="LC1" |colData(sce)$Identity=="LC2" |colData(sce)$Identity=="LP")]

#mapping cell types onto sce
Type <- c("LC","LC","LP")
Identity <- c("LC1","LC2","LP")
colData(sce.sub)$Type <- mapvalues(colData(sce.sub)$Identity, Identity, Type)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sce.sub <- runDiffusionMap(sce.sub, dimred="corrected")
dm <- reducedDim(sce.sub, "DiffusionMap")
write.csv(dm, file= "../../data/Fig_4/DiffusionCoordinates_dimredCorrected.csv")

#Supplementary_Fig_15
plotReducedDim(sce.sub, dimred="DiffusionMap", colour_by="Identity") + scale_color_manual(values=setNames(colours.classes, classes))
ggsave(filename="../../data/Supp_Fig_15/Supp_Fig_15_i.pdf",width=8,height=7)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(TSCAN)
pseudo.out <- quickPseudotime(sce.sub, sce.sub$Identity, use.dimred="corrected", outgroup=TRUE)
common.pseudo <- rowMeans(pseudo.out$ordering, na.rm=TRUE)

plotReducedDim(sce.sub, colour_by=I(common.pseudo), dimred="DiffusionMap")
ggsave(filename="../../data/Supp_Fig_15/Supp_Fig_15_ii.pdf",width=8,height=7)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

saveRDS(sce.sub,file="../../data/sce_Luminal_DiffMap.rds")
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sessionInfo()

