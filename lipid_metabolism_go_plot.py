# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# ===================== Lipid Metabolism GO Enrichment Bubble Plot =====================
file_path = r"D:\ZHB\HBD-2026\脂代谢GO\_FINAL_GO.csv"
df = pd.read_csv(file_path)

df["LogP_abs"] = df["LogP"].abs()

go_list = [
    ("GO Biological Processes", "Biological Process"),
    ("GO Cellular Components", "Cellular Component"),
    ("GO Molecular Functions", "Molecular Function")
]

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8

for cat_name, full_name in go_list:
    sub_df = df[df["Category"] == cat_name].copy()
    if sub_df.empty:
        continue

    sub_df = sub_df.sort_values("LogP_abs", ascending=False).head(10)

    plt.figure(figsize=(7, 5))

    sns.scatterplot(
        data=sub_df,
        y="Description",
        x="Enrichment",
        size="#GeneInGOAndHitList",
        hue="LogP_abs",
        palette="viridis",
        sizes=(300, 1200),
        edgecolor="black",
        linewidth=0.8
    )

    plt.title(f"GO Enrichment - {full_name}", fontsize=11, fontweight='bold')
    plt.xlabel("Enrichment Factor", fontsize=9)
    plt.ylabel("")
    plt.grid(alpha=0.2, linestyle='--')
    plt.tight_layout()

    save_path = rf"D:\ZHB\HBD-2026\{full_name.replace(' ','_')}.svg"
    plt.savefig(save_path, format="svg", dpi=300, bbox_inches="tight")
    print(f"✅ Saved: {save_path}")
    plt.show()