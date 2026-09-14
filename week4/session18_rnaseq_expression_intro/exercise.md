# Session 18 — Exercise: Full differential expression pass

Use `data/arabidopsis_drought_rnaseq_counts.csv` (20 genes, 3 control + 3 drought replicates).

## Tasks
1. Run the code-along. Report the top 3 upregulated and top 3 downregulated genes by log2FC.
2. Compute the **coefficient of variation** (std/mean) across replicates for each gene's control samples,
   as a crude "how noisy is this gene's measurement" check. Are the always-on "housekeeping" genes
   (ACTIN2, GAPDH, UBQ10, TUB4) less variable than the stress-response genes? (They should be — that's
   exactly why they're chosen as normalization/reference genes in real qPCR experiments!)
3. Using the `|log2FC| > 1` threshold from the code-along, list all "differentially expressed" genes and
   cross-reference each with the GO annotation dictionary from Session 17 — do they mostly fall into
   "response to water deprivation" or "photosynthesis" categories, matching what we'd expect?
4. **Discussion:** the code-along uses a fold-change-only threshold, with no statistical test. A gene could
   have a huge log2FC just from one noisy outlier replicate. What additional check (using the 3 replicates
   you have) could give you more confidence a gene's change is real, short of a full statistical test?
5. **Stretch:** build a simple "confidence" score = |log2FC| / (variability across drought replicates), and
   re-rank genes by this instead of raw log2FC alone. Does the gene ranking change much?

## Answer key
```python
import pandas as pd, numpy as np

counts = pd.read_csv("../../data/arabidopsis_drought_rnaseq_counts.csv", index_col="gene_id")
cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
control_cols = ["Control_rep1","Control_rep2","Control_rep3"]

cv = cpm[control_cols].std(axis=1) / cpm[control_cols].mean(axis=1)
print(cv.sort_values())
# Housekeeping genes (ACTIN2, GAPDH, UBQ10, TUB4) should show noticeably lower CV
# than genes whose true biological expression is changing a lot.
```
