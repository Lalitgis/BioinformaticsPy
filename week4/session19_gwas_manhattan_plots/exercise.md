# Session 19 — Exercise: Fine-map a GWAS peak

Use `data/rice_plant_height_gwas.csv`.

## Tasks
1. Run the code-along and view `manhattan_plot.png`. Identify visually which chromosome has the tallest peak.
2. Report the top SNP's exact position and p-value, and define a "candidate region" as ± 100,000 bp around
   it. List all SNPs that fall in that window and are also individually significant (p < 5e-8).
3. Compute how many chromosomes show **zero** significant SNPs versus how many show at least one — this is
   what a real GWAS with a true single strong locus usually looks like (one clear peak, flat elsewhere).
4. **Discussion:** if a breeder wanted to develop a CAPS marker (Session 13!) to track this height-associated
   region during breeding, what would they need in addition to the GWAS peak's position? (Hint: they need
   the actual DNA sequence at that position, and ideally a SNP that also happens to fall in/create a
   restriction site.)
5. **Stretch:** re-run the whole analysis using a stricter significance threshold (`1e-10` instead of
   `5e-8`). How does the count of "significant" SNPs change? Discuss the tradeoff between a stricter
   threshold (fewer false positives) and a looser one (fewer missed true associations, i.e. false negatives).

## Answer key
```python
import pandas as pd

gwas = pd.read_csv("../../data/rice_plant_height_gwas.csv")
top_hit = gwas.loc[gwas["p_value"].idxmin()]
window = gwas[(gwas["chrom"] == top_hit["chrom"]) &
              (gwas["pos"].between(top_hit["pos"] - 100000, top_hit["pos"] + 100000))]
sig_in_window = window[window["p_value"] < 5e-8]
print(f"{len(sig_in_window)} significant SNPs within 100kb of the top hit")
print(sig_in_window[["snp_id", "pos", "p_value"]].sort_values("pos"))
```
