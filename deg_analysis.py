import os
import pandas as pd
import numpy as np
from scipy import stats
import glob
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

data_dir = r"D:\ZHB\HBD-2026\GSE98969\GSE98969_RAW (2)"
design_file = os.path.join(data_dir, "design.txt")

design_df = pd.read_csv(design_file, sep='\t')

def map_group(x):
    x_str = str(x).strip()
    if '5XFAD' in x_str:
        return 'AD'
    elif 'C57BL/6' in x_str or 'WT' in x_str:
        return 'Control'
    else:
        return 'Other'

design_df['Group'] = design_df['Mouse_ID'].apply(map_group)
design_df = design_df[design_df['Group'].isin(['Control', 'AD'])]
cell_groups = design_df.set_index('Well_ID')['Group']

print("正在读取数据...")
all_txt = glob.glob(os.path.join(data_dir, "*.txt"))
expr_files = [f for f in all_txt if "design.txt" not in f]

df_list = []
for f in tqdm(expr_files, desc="读取文件"):
    try:
        tmp = pd.read_csv(f, sep='\t', index_col=0)
        df_list.append(tmp)
    except:
        continue

expr_df = pd.concat(df_list, axis=1)
common_cells = expr_df.columns.intersection(cell_groups.index)
expr_df = expr_df[common_cells]
cell_groups = cell_groups[common_cells]

print(f"✅ 读取完成：基因 {expr_df.shape[0]}，细胞 {expr_df.shape[1]}")

ctrl_cells = cell_groups[cell_groups == "Control"].index
ad_cells = cell_groups[cell_groups == "AD"].index

res = []
for gene in tqdm(expr_df.index, desc="差异分析"):
    x = expr_df.loc[gene, ctrl_cells].astype(float).dropna()
    y = expr_df.loc[gene, ad_cells].astype(float).dropna()
    if len(x) < 3 or len(y) < 3:
        continue

    log2fc = np.log2(y.mean() + 1e-8) - np.log2(x.mean() + 1e-8)
    _, p = stats.mannwhitneyu(x, y, alternative='two-sided')

    res.append({
        "Gene": gene,
        "Control_mean": x.mean(),
        "AD_mean": y.mean(),
        "log2FC": log2fc,
        "p_value": p
    })

res_df = pd.DataFrame(res)
res_df = res_df.sort_values("p_value")
res_df["FDR"] = multipletests(res_df["p_value"], method="fdr_bh")[1]

sig = res_df[(res_df["FDR"] < 0.05) & (res_df["log2FC"].abs() > 0.25)]

print(f"\n显著差异基因：{len(sig)}")
res_df.to_excel(os.path.join(data_dir, "All_DEGs.xlsx"), index=False)
sig.to_excel(os.path.join(data_dir, "Significant_DEGs.xlsx"), index=False)

print("\n🎉 全部完成！")