# LOOCV.py —— 适用于小数表达矩阵 · 最终无错版
import pandas as pd
import numpy as np
import os
from scipy import stats
import warnings

warnings.filterwarnings('ignore')

# ===================== 路径 =====================
excel_file = r"D:\ZHB\HBD-2026\LOOCV\LOOCV counts.xlsx"
out_csv = r"D:\ZHB\HBD-2026\LOOCV\LOOCV_reproducibility.csv"


# ===================== 工具函数 =====================
def parse_condition(col):
    return col.split('_')[0]


# 差异分析：t检验 · 完美兼容小数表达矩阵
def deg_analysis(expr_df, meta, lfc_thr=0.5, p_thr=0.05):
    genes = expr_df.index
    res_list = []

    # 提取分组
    samples_apoe3 = meta[meta["condition"] == "APOE3"]["sample"].values
    samples_apoe4 = meta[meta["condition"] == "APOE4"]["sample"].values

    for gene in genes:
        # 提取该基因在两组中的表达值
        expr3 = expr_df.loc[gene, samples_apoe3].astype(float).values
        expr4 = expr_df.loc[gene, samples_apoe4].astype(float).values

        # 计算倍数变化 & t检验
        mean3 = np.mean(expr3)
        mean4 = np.mean(expr4)
        lfc = np.log2((mean4 + 1e-6) / (mean3 + 1e-6))
        t_stat, p_val = stats.ttest_ind(expr4, expr3, equal_var=False)

        # 判断是否差异基因
        is_deg = (p_val < p_thr) and (abs(lfc) > lfc_thr)
        res_list.append({"gene": gene, "log2FC": lfc, "pvalue": p_val, "is_DEG": is_deg})

    return pd.DataFrame(res_list).set_index("gene")


# ===================== 读取数据 =====================
counts = pd.read_excel(excel_file, index_col=0)
# 去除重复基因（防止报错）
counts = counts[~counts.index.duplicated(keep="first")]

samples = counts.columns.tolist()
meta = pd.DataFrame({
    "sample": samples,
    "condition": [parse_condition(s) for s in samples]
})

n = len(samples)
deg_count = {g: 0 for g in counts.index}

# ===================== LOOCV 核心循环 =====================
print("=" * 65)
print("        LOOCV 留一法交叉验证（小数表达矩阵专用）")
print("=" * 65)

for i, leave_out in enumerate(samples):
    print(f"\n【第 {i + 1}/{n} 轮】剔除样本: {leave_out}")

    # 保留其余样本
    keep_samples = [s for s in samples if s != leave_out]
    expr_sub = counts[keep_samples]
    meta_sub = meta[meta["sample"].isin(keep_samples)].copy()

    # 必须保留两组才能做差异分析
    if len(meta_sub["condition"].unique()) < 2:
        print("  → 跳过：只剩一组样本")
        continue

    try:
        res = deg_analysis(expr_sub, meta_sub)
        degs = res[res["is_DEG"]].index.tolist()
        print(f"  → 成功！识别差异基因: {len(degs)} 个")

        for g in degs:
            if g in deg_count:
                deg_count[g] += 1

    except Exception as e:
        print(f"  → 错误: {str(e)[:80]}")

# ===================== 计算重现率 =====================
repro_rate = {g: round((cnt / n) * 100, 2) for g, cnt in deg_count.items()}

# ===================== 输出关键基因 =====================
key_genes = ["NRF1", "PPARA", "UCP3", "PDK4", "SIRT1", "SIRT3", "FOXO1", "FOXO3", "SOD2"]

print("\n" + "=" * 65)
print("                   关键基因重现率")
print("=" * 65)
for g in key_genes:
    rate = repro_rate.get(g, 0.0)
    print(f"  {g:10s} →  {rate:>6.1f} %")

# ===================== 保存结果 =====================
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
df_result = pd.DataFrame(list(repro_rate.items()), columns=["gene", "reproducibility"])
df_result.to_csv(out_csv, index=False)

print("\n✅ 全部运行完成！结果已保存到：")
print(f"   {out_csv}")