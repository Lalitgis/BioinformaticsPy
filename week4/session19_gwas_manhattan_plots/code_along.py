"""
Session 19 — GWAS concepts & Manhattan plots
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

gwas = pd.read_csv("../../data/rice_plant_height_gwas.csv")
gwas["p_value"] = gwas["p_value"].astype(float)
gwas["neg_log10_p"] = -np.log10(gwas["p_value"])

print(f"Loaded {len(gwas)} SNPs across {gwas['chrom'].nunique()} chromosomes")
print(gwas.sort_values("p_value").head(5))

# --- Order for plotting: by chromosome, then position ---
gwas = gwas.sort_values(["chrom", "pos"]).reset_index(drop=True)
gwas["plot_pos"] = range(len(gwas))

chrom_centers = gwas.groupby("chrom")["plot_pos"].median()

plt.figure(figsize=(13, 4.5))
colors = ["#4a90d9", "#c98b3a"]
for i, chrom in enumerate(sorted(gwas["chrom"].unique())):
    subset = gwas[gwas["chrom"] == chrom]
    plt.scatter(subset["plot_pos"], subset["neg_log10_p"], s=6, color=colors[i % 2])

sig_threshold = -np.log10(5e-8)
plt.axhline(sig_threshold, color="red", linestyle="--", linewidth=1, label="genome-wide significance (5e-8)")
plt.xticks(chrom_centers.values, [f"Chr{c}" for c in chrom_centers.index], rotation=90, fontsize=7)
plt.xlabel("Chromosome")
plt.ylabel("-log10(p-value)")
plt.title("Rice plant-height GWAS")
plt.legend()
plt.tight_layout()
plt.savefig("manhattan_plot.png", dpi=150)
print("\nSaved manhattan_plot.png")

# --- Report the top hit ---
top_hit = gwas.loc[gwas["p_value"].idxmin()]
print(f"\nTop associated SNP: {top_hit['snp_id']} on chromosome {top_hit['chrom']}, "
      f"position {top_hit['pos']}, p-value = {top_hit['p_value']:.2e}")

# --- Count SNPs passing genome-wide significance ---
sig_snps = gwas[gwas["p_value"] < 5e-8]
print(f"{len(sig_snps)} SNPs pass the genome-wide significance threshold")
print(f"They span chromosome(s): {sorted(sig_snps['chrom'].unique())}, "
      f"positions {sig_snps['pos'].min()}-{sig_snps['pos'].max()}")
