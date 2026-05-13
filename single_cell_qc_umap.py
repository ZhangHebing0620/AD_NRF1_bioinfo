# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

# ===================== Single-cell QC & UMAP Analysis =====================
plt.rcParams['savefig.format'] = 'svg'
plt.rcParams['figure.dpi'] = 300

file_path = r"D:\ZHB\HBD-2026\GSE138852\GSE138852_counts.csv"
print("Reading data...")
adata = sc.AnnData(pd.read_csv(file_path, index_col=0).T)

# QC metrics
adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1)
adata.obs['total_counts'] = adata.X.sum(axis=1)

# Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalization & scaling
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata.raw = adata

# Dimensionality reduction
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# QC violin plot
print("Generating QC violin plot...")
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts'],
    jitter=0.4,
    multi_panel=True,
    save="_QC.svg",
    linewidth=0,
    alpha=0.3,
    color="#6a9ace"
)

# Visualization
sc.pl.umap(adata, title="UMAP Clustering", save="_UMAP.svg")
sc.pl.umap(adata, color="NRF1", cmap="Blues", title="NRF1 Expression", save="_NRF1.svg")
sc.pl.umap(adata, color=["SNAP25", "GFAP", "PECAM1"], ncols=3, save="_markers.svg")

print("✅ All figures saved in 'figures' folder")