#!/bin/bash

# Make sure the right conda env is activated
if ! command -v pyscenic &> /dev/null
then
        echo "PYSCENIC could not be found, activate conda env"
	    exit
	fi

LOOM_FILE="../data/Loom_LumAllsub_subset1k.loom"
TF_FILE="../data/forPySCENIC/hs_hgnc_tfs.txt"
REG="../data/SCENIC_reg_subset1k.csv"
OUT="../data/SCENIC_OUT_subset1k.loom"
LOGS="../logs/"

if test -f "$LOGS/AUC_subset1k.out"; then
    rm $LOGS/AUC_subset1k.out
    rm $LOGS/AUC_subset1k.err
fi

bsub -o $LOGS/AUC_subset1k.out -e $LOGS/AUC_subset1k.err -M 120000 -n 24 pyscenic aucell $LOOM_FILE $REG --out $OUT --num_workers 20
