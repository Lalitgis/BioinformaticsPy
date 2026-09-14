"""
Session 18 — RNA-seq expression analysis basics
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

counts = pd.read_csv("../../data/arabidopsis_drought_rnaseq_counts.csv", index_col="gene_id")
print("Raw counts:")
print(counts)

cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6

control_cols = ["Control_rep1", "Control_rep2", "Control_rep3"]
drought_cols = ["Drought_rep1", "Drought_rep2", "Drought_rep3"]

mean_control = cpm[control_cols].mean(axis=1)
mean_drought = cpm[drought_cols].mean(axis=1)
log2fc = np.log2((mean_drought + 1) / (mean_control + 1))

results = pd.DataFrame({
    "mean_control_cpm": mean_control.round(1),
    "mean_drought_cpm": mean_drought.round(1),
    "log2FC": log2fc.round(2),
})
results_sorted = results.sort_values("log2FC", ascending=False)
print("\nGenes ranked by log2 fold change (Drought vs Control):")
print(results_sorted)

# --- Simple "differentially expressed" call: |log2FC| > 1 (i.e., >= 2-fold change) ---
de_genes = results_sorted[results_sorted["log2FC"].abs() > 1]
print(f"\n{len(de_genes)} of {len(results)} genes pass the |log2FC| > 1 threshold:")
print(de_genes)

# --- Bar plot of top up/down genes ---
top_up = results_sorted.head(6)
top_down = results_sorted.tail(6)
combined = pd.concat([top_up, top_down])

plt.figure(figsize=(8, 5))
colors = ["#4a9e3f" if v > 0 else "#c9483a" for v in combined["log2FC"]]
plt.barh(combined.index, combined["log2FC"], color=colors)
plt.xlabel("log2 fold change (Drought vs Control)")
plt.title("Top up/down-regulated genes under drought stress")
plt.axvline(0, color="black", linewidth=0.8)
plt.tight_layout()
plt.savefig("top_drought_genes.png", dpi=150)
print("\nSaved top_drought_genes.png")
