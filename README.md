# LactatingMammaryCells

## About
This repository contains all scripts to reproduce the results from the paper "Transcriptional changes in the mammary gland during lactation revealed by single cell sequencing of cells from human milk".

## Overview
![](MilkCells.png)

## Abstract
Findings from epidemiological studies have shown that breast cancer risk is influenced by parity in an age-dependent manner. However, human mammary tissue remodelling that takes place during pregnancy and lactation remains poorly understood due to the challenge of acquiring samples. We report here single-cell transcriptomic analysis of 110,744 viable breast cells isolated from human milk or non-lactating breast tissue, isolated from nine and seven donors, respectively. We found that human milk largely contains epithelial cells belonging to the luminal lineage and a repertoire of immune cells. Further transcriptomic analysis of the milk cells identified two distinct secretory cell types that shared similarities with luminal progenitors, but no populations comparable to hormone-responsive cells. Taken together, our data offers a comprehensive reference map and a window on the cellular dynamics that occur during human lactation and provides further insights on the interplay between pregnancy, lactation and breast cancer.

## Repository
The repository is structured as follows:

- All scripts are contained in the [src](src/) folder
- Initial analysis was conducted on the [Individual_batches](src/1_Individual_batches) before down stream analysis
- To run downstream analysis run the scripts (in order) in the [Preparing_merged_data](src/2_Preparing_merged_data) to generate the final sce file
- To replicate any of the analysis in the [Figures](src/3a_Figures) or [Supplementary_Figures](src/3b_Supplementary_Figs) follow the scripts indicated. Please note that in some special cases some scripts must be run prior, however these are indicated at the top of the script.

## Sample overview
Overall samples from non-lactation associated mammary cells (NMC) were isolated from resting breast (RB) tissues and lactation associated mammary cells (LMC) consisted of human milk cells (HMC). Samples were sequenced in 3 separate batches. Batch 1 contains samples LMC1-4 and NMC1-4; Batch 2 contains samples LMC2B, LMC5-8 and NMC5-7; Batch 3 contains samples LMC9 and NMC1B.

## Links
- Interactive [website](http://bioinf.stemcells.cam.ac.uk:3838/khaled_wUFt1bHfmC/twigger/) to explore data online
- Data can also be accessed at Array Express: 
Batch 1: ([E-MTAB-9841](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9841))
Batch 2: ([E-MTAB-10855](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10855))
Batch 3: ([E-MTAB-10885](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10885))
- [WTK Lab](https://www.phar.cam.ac.uk/research/Khaled)

## Questions
For questions contact "ajt215@cam.ac.uk"