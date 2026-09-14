# Session 17 — Exercise: Manual GO enrichment on the drought gene panel

Use the `go_annotation` dictionary from the code-along.

## Tasks
1. Run the code-along and confirm "response to water deprivation" is strongly enriched in the
   drought-upregulated gene list.
2. Build a "drought-downregulated" gene list yourself (photosynthesis genes should be it — check
   Session 18's dataset design if you want to peek ahead) and check enrichment of "photosynthesis, light
   reaction" in that list instead.
3. Add 2 more GO-style terms of your own choosing to the `go_annotation` dictionary for any 2 genes not yet
   fully annotated (e.g., give `RBCS1A` an additional Cellular Component of `"chloroplast"`), then extend
   `enrichment_check()` to work on any namespace (BP, MF, or CC), not just BP.
4. **Discussion:** our enrichment fold numbers here are simple ratios, with no statistical test. What could
   go wrong if you reported "10x enrichment!" for a GO term that only had 1 gene in a background of 20 —
   compare that to 1 gene in a background of 20,000. Why does list/background *size* matter for how much
   you should trust an enrichment result?
5. **Stretch:** look up (or design based on real biological reasoning) what real GO Biological Process term
   `ABI5` and `NCED3` would likely share, beyond "response to abscisic acid" — both are central players in
   ABA (abscisic acid) hormone signaling during drought stress.

## Answer key
```python
drought_down = ["RBCS1A", "RBCS2B", "CAB1", "PSBA", "PIP2_1", "EXPA1"]

def enrichment_check(gene_list, term, all_genes, namespace="BP"):
    hits_in_list = sum(1 for g in gene_list if go_annotation.get(g, {}).get(namespace) == term)
    hits_in_bg = sum(1 for g in all_genes if go_annotation.get(g, {}).get(namespace) == term)
    return hits_in_list / len(gene_list), hits_in_bg / len(all_genes)

frac_l, frac_b = enrichment_check(drought_down, "photosynthesis, light reaction", all_genes)
print("In list:", frac_l, "in background:", frac_b)
```
Discussion: with tiny sample sizes (1 gene out of 20), a single gene shifting categories can create a huge
apparent "fold enrichment" purely by chance — this is exactly why real tools use a statistical test
(Fisher's exact test / hypergeometric test) that accounts for sample size, not a raw ratio, and why you
should be skeptical of enrichment claims based on very short gene lists.
