# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

plt.rcParams["savefig.format"] = "svg"
plt.rcParams["savefig.dpi"] = 300

file_path = r"D:\ZHB\HBD-2026\GSE138852\GSE138852_counts.csv"

print("正在读取数据...")
counts = pd.read_csv(file_path, index_col=0)
adata = sc.AnnData(counts.T)
print("数据维度：", adata.shape)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(adata, title="UMAP distribution", save="_UMAP.svg", show=False)
sc.pl.umap(adata, color="NRF1", cmap="Blues", title="NRF1 expression", save="_NRF1_UMAP.svg", show=False)

print("✅ UMAP 分析完成！")