# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# GO Enrichment Bubble Plot for NRF1-related genes
file_path = r"D:\ZHB\HBD-2026\GO-30-NRF1\Enrichment_GO\_FINAL_GO.csv"
df = pd.read_csv(file_path)

print("All Categories in file:")
print(df["Category"].unique())

df["LogP_abs"] = df["LogP"].abs()

# GO categories
go_list = [
    ("GO Biological Processes", "BP"),
    ("GO Cellular Components", "CC"),
    ("GO Molecular Functions", "MF")
]

for cat_name, label in go_list:
    sub_df = df[df["Category"] == cat_name].copy()
    if sub_df.empty:
        print(f"⚠️ {label} No data, skip")
        continue

    sub_df = sub_df.sort_values("LogP_abs", ascending=False).head(10)

    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=sub_df,
        y="Description",
        x="Enrichment",
        size="#GeneInGOAndHitList",
        hue="LogP_abs",
        palette="RdYlBu_r",
        sizes=(200, 1000),
        edgecolor="black"
    )
    plt.title(f"GO Enrichment - {label}", fontsize=14, fontweight='bold')
    plt.tight_layout()

    save_path = rf"D:\ZHB\HBD-2026\GO_{label}.svg"
    plt.savefig(save_path, format="svg", bbox_inches="tight")
    print(f"✅ Saved: {save_path}")
    plt.show()