library(SingleCellExperiment)
library(reticulate)
sce <- readRDS("../data/sce_LumAllsub.rds")
rD <- rowData(sce)
rD <- rD[,c("Symbol","hvg")]
#Adata doesn't like the rotation matrix in the rD
sce <- SingleCellExperiment(assays=SimpleList("counts"=counts(sce),
					      "logcounts"=logcounts(sce)),
			    colData=colData(sce),
			    rowData=rD)
use_condaenv('base')
sceasy::convertFormat(sce, from="sce", to="anndata",
                       outFile='../data/Adata_LumAllsub.h5ad')
