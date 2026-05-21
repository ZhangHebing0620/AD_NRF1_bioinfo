import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==============================================
# 1. 你的数据集信息（样本量 + logFC阈值）
# ==============================================
datasets = {
    "D1": {"n": 863,   "logFC_cutoff": 1.0},
    "D2": {"n": 5234,  "logFC_cutoff": 1.0},
    "D3": {"n": 1303,  "logFC_cutoff": 0.2},
    "D4": {"n": 2844,  "logFC_cutoff": 0.2}
}

# ==============================================
# 2. 自动计算权重（样本量 + 严格度）
# ==============================================
def calculate_weight(n_sample, logFC_cutoff):
    w_sample = np.sqrt(n_sample)
    w_strict = 1.5 if logFC_cutoff >= 0.8 else 0.5
    final_w = w_sample * w_strict
    return round(final_w, 1)

w1 = calculate_weight(datasets["D1"]["n"], datasets["D1"]["logFC_cutoff"])
w2 = calculate_weight(datasets["D2"]["n"], datasets["D2"]["logFC_cutoff"])
w3 = calculate_weight(datasets["D3"]["n"], datasets["D3"]["logFC_cutoff"])
w4 = calculate_weight(datasets["D4"]["n"], datasets["D4"]["logFC_cutoff"])

print("✅ 四个数据集最终权重：")
print(f"D1: {w1}")
print(f"D2: {w2}")
print(f"D3: {w3}")
print(f"D4: {w4}\n")

# ==============================================
# 3. 你的文件绝对路径（我已经帮你写好）
# ==============================================
base_path = r"D:\ZHB\HBD-2026\加权投票机制"

deg1 = set(pd.read_csv(f"{base_path}\\deg1.csv", header=0)["gene"])
deg2 = set(pd.read_csv(f"{base_path}\\deg2.csv", header=0)["gene"])
deg3 = set(pd.read_csv(f"{base_path}\\deg3.csv", header=0)["gene"])
deg4 = set(pd.read_csv(f"{base_path}\\deg4.csv", header=0)["gene"])

# ==============================================
# 4. 加权投票打分
# ==============================================
all_genes = sorted(deg1 | deg2 | deg3 | deg4)
results = []

for gene in all_genes:
    score = 0
    if gene in deg1: score += w1
    if gene in deg2: score += w2
    if gene in deg3: score += w3
    if gene in deg4: score += w4

    results.append({
        "gene": gene,
        "in_D1": 1 if gene in deg1 else 0,
        "in_D2": 1 if gene in deg2 else 0,
        "in_D3": 1 if gene in deg3 else 0,
        "in_D4": 1 if gene in deg4 else 0,
        "weighted_score": score
    })

df = pd.DataFrame(results)
df = df.sort_values("weighted_score", ascending=False)

# 保存结果到同一个文件夹
df.to_csv(f"{base_path}\\weighted_voting_result.csv", index=False)
print(f"✅ 结果已保存到：{base_path}\\weighted_voting_result.csv")

# ==============================================
# 5. 画图：得分分布图
# ==============================================
plt.figure(figsize=(7, 4))
plt.hist(df["weighted_score"], bins=30, color="#4a90e2", edgecolor="black", alpha=0.8)
plt.xlabel("Weighted Consensus Score")
plt.ylabel("Number of Genes")
plt.title("Weighted Voting Score Distribution")
plt.tight_layout()
plt.savefig(f"{base_path}\\weighted_score_distribution.svg", dpi=300)
plt.show()

print("✅ 分布图已保存：weighted_score_distribution.svg")