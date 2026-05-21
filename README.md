# AD_NRF1_bioinfo
Reproducible code for bioinformatics analysis of NRF1 in Alzheimer's disease

## Overview
This repository contains Python scripts for all analyses in this study:
1.  Single-cell RNA-seq QC, UMAP visualization, and marker gene plotting
2.  Differential gene expression analysis
3.  GO enrichment bubble plot for NRF1-related genes
4.  GO enrichment bubble plot for lipid metabolism-related genes
5.  Clinical correlation analysis of NRF1 expression with Braak/CDR/CERAD scores
6.  MCODE module identification and K-Core parameter sensitivity analysis
7.  Weighted voting algorithm for multi-dataset consensus gene selection
8.  Leave-One-Out Cross-Validation (LOOCV) for model stability

## Usage
Run each script with Python 3:
```bash
python single_cell_qc_umap.py
python deg_analysis.py
python go_enrichment_plot.py
python lipid_metabolism_go_plot.py
python clinical_correlation.py
python mcode_stability.py
python weighted_voting.py
python LOOCV.py
