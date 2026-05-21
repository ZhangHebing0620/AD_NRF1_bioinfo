# -*- coding: utf-8 -*-
"""
Graph Propagation Analysis for PPI Network
Integrate weighted voting scores with protein-protein interaction network
For AD_NRF1_bioinfo project
"""
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8

# ===================== 1. Read STRING PPI data =====================
def read_ppi_network(file_path, score_threshold=0.7):
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
    if first_line.startswith('#'):
        columns = first_line.lstrip('#').split('\t')
    else:
        columns = None

    df = pd.read_csv(file_path, sep='\t', skiprows=1, names=columns)
    score_col = 'combined_score'
    if score_col not in df.columns:
        raise ValueError(f"Column {score_col} not found in PPI file")

    df_filtered = df[df[score_col] >= score_threshold].copy()
    if len(df_filtered) == 0:
        df_filtered = df.copy()
    G = nx.from_pandas_edgelist(df_filtered, 'node1', 'node2')
    return G

# ===================== 2. Read weighted voting result =====================
def read_weighted_voting(file_path):
    df = pd.read_csv(file_path)
    if 'gene' not in df.columns:
        if 'Gene' in df.columns:
            df.rename(columns={'Gene': 'gene'}, inplace=True)
        elif 'Gene name' in df.columns:
            df.rename(columns={'Gene name': 'gene'}, inplace=True)
    df['gene'] = df['gene'].str.upper()
    df['rank'] = df['weighted_score'].rank(ascending=False, method='dense')
    df['f0'] = 1.0 / df['rank']
    return df

# ===================== 3. Graph propagation algorithm =====================
def graph_propagation(G, f0, alpha=0.5, max_iter=100, tol=1e-6):
    nodes = list(G.nodes)
    f = f0.copy()
    for iteration in range(max_iter):
        f_new = {}
        for node in nodes:
            neighbors = list(G.neighbors(node))
            if neighbors:
                neighbor_mean = np.mean([f.get(nei, 0) for nei in neighbors])
            else:
                neighbor_mean = f.get(node, 0)
            f_new[node] = alpha * neighbor_mean + (1 - alpha) * f0.get(node, 0)
        diff = max(abs(f_new[n] - f.get(n, 0)) for n in nodes)
        f = f_new
        if diff < tol:
            print(f"Converged at iteration {iteration+1}, diff={diff:.2e}")
            break
    return f

# ===================== 4. Main function =====================
if __name__ == "__main__":
    # You can modify paths when running
    ppi_file = "string_interactions_short-30-2.tsv"
    vote_file = "weighted_voting_result.csv"
    output_dir = "./"

    # Load network and genes
    G = read_ppi_network(ppi_file)
    df_vote = read_weighted_voting(vote_file)

    # Match genes in network
    G = nx.relabel_nodes(G, lambda x: str(x).upper())
    genes_in_graph = set(G.nodes)
    df_vote = df_vote[df_vote['gene'].isin(genes_in_graph)].copy()
    f0_dict = dict(zip(df_vote['gene'], df_vote['f0']))

    # Run propagation
    f_prop = graph_propagation(G, f0_dict, alpha=0.5)

    # Rank results
    df_prop = pd.DataFrame(list(f_prop.items()), columns=['gene', 'prop_score'])
    df_prop = df_prop.sort_values('prop_score', ascending=False)
    df_prop['prop_rank'] = np.arange(1, len(df_prop)+1)

    # Merge and compare
    df_final = df_vote.merge(df_prop, on='gene', how='left')
    df_final['rank_change'] = df_final['prop_rank'] - df_final['rank']

    # Output NRF1 and module genes
    print("\n=== NRF1 Ranking ===")
    nrf1 = df_final[df_final['gene'] == 'NRF1']
    if not nrf1.empty:
        print(f"Original rank: {nrf1['rank'].iloc[0]}, Propagated rank: {nrf1['prop_rank'].iloc[0]}")

    # Lipid metabolism module
    module_genes = ['PPARA', 'UCP3', 'PDK4', 'SIRT1', 'SIRT3', 'FOXO1', 'FOXO3', 'SOD2']
    module_df = df_final[df_final['gene'].isin(module_genes)]
    print("\n=== Module Genes Ranking Change ===")
    print(module_df[['gene', 'rank', 'prop_rank', 'rank_change']])

    # Save results
    df_final.to_csv(f"{output_dir}/graph_propagation_results.csv", index=False)
    module_df.to_csv(f"{output_dir}/module_genes_propagation.csv", index=False)

    print("\n✅ Graph propagation analysis completed!")