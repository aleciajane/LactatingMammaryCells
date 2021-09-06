#!/bin/bash

# Make sure the right conda env is activated
if ! command -v pyscenic &> /dev/null
then
        echo "PYSCENIC could not be found, activate conda env"
	    exit
	fi

LOOM_FILE="../data/Loom_LumAllsub_subset1k.loom"
TF_FILE="../data/forPySCENIC/hs_hgnc_tfs.txt"
ADJ="../data/SCENIC_adj_subset1k.csv"
FEATHER1="../data/forPySCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
FEATHER2="../data/forPySCENIC/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
ANNOS="../data/forPySCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
OUT="../data/SCENIC_reg_subset1k.csv"
LOGS="../logs/"

if test -f "$LOGS/REG_subset1k.out"; then
    rm $LOGS/REG_subset1k.out
    rm $LOGS/REG_subset1k.err
fi

bsub -o $LOGS/REG_subset1k.out -e $LOGS/REG_subset1k.err -M 80000 -n 24 pyscenic ctx $ADJ $FEATHER1 $FEATHER2 --annotations_fname $ANNOS --expression_mtx_fname $LOOM_FILE --mode "dask_multiprocessing" --output $OUT --num_workers 20 --mask_dropouts
