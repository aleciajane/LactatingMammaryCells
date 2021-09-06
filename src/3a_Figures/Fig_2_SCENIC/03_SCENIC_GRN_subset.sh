#!/bin/bash

# Make sure the right conda env is activated
if ! command -v pyscenic &> /dev/null
then
        echo "PYSCENIC could not be found, activate conda env"
	    exit
	fi

LOOM_FILE="../data/Loom_LumAllsub_subset1k.loom"
TF_FILE="../data/forPySCENIC/hs_hgnc_tfs.txt"
OUT="../data/SCENIC_adj_subset1k.csv"
LOGS="../logs/"

if test -f "$LOGS/GRN_subset1k.out"; then
    rm $LOGS/GRN_subset1k.out
    rm $LOGS/GRN_subset1k.err
fi

bsub -o $LOGS/GRN_subset1k.out -e $LOGS/GRN_subset1k.err -M 80000 -n 24 pyscenic grn $LOOM_FILE $TF_FILE -o $OUT --num_workers 20
