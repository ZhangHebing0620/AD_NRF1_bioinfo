def get_genes_from_sif(sif_path):
    genes = set()
    with open(sif_path, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                genes.add(parts[0])
                genes.add(parts[2])
    return genes

def jaccard_index(set_a, set_b):
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union != 0 else 0

if __name__ == "__main__":
    cluster2 = get_genes_from_sif("MCODE_default_Cluster2.sif")
    cluster3 = get_genes_from_sif("MCODE_default_Cluster3.sif")
    jacc = jaccard_index(cluster2, cluster3)
    print(f"Jaccard = {jacc:.3f}")
    print(f"Cluster2 genes: {len(cluster2)}")
    print(f"Cluster3 genes: {len(cluster3)}")