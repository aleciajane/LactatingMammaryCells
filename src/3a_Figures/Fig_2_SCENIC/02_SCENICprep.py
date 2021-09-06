import pandas as pd
import scanpy as sc
import loompy as lp
import os
import numpy as np

sce = sc.read_h5ad("../data/Adata_LumAllsub.h5ad")
sc.pp.filter_genes(sce, min_cells=3)

# create basic row and column attributes for the loom file:
f_loom_path_scenic="../data/Loom_LumAllsub.loom"
row_attrs = {
    "Gene": np.array(sce.var_names) ,
}
col_attrs = {
    "CellID": np.array(sce.obs_names) ,
    "nGene": np.array( np.sum(sce.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(sce.X.transpose() , axis=0)).flatten() ,
}
lp.create(f_loom_path_scenic, sce.X.transpose(), row_attrs, col_attrs)

rnd_cells = np.random.choice(np.arange(1,sce.n_obs),size=20000,replace=False)
sce = sce[rnd_cells]
# create basic row and column attributes for the loom file:
f_loom_path_scenic="../data/Loom_LumAllsub_subset1k.loom"
row_attrs = {
    "Gene": np.array(sce.var_names) ,
}
col_attrs = {
    "CellID": np.array(sce.obs_names) ,
    "nGene": np.array( np.sum(sce.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(sce.X.transpose() , axis=0)).flatten() ,
}
lp.create(f_loom_path_scenic, sce.X.transpose(), row_attrs, col_attrs)
