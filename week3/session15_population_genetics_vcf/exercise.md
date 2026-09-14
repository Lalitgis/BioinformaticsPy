# Session 15 — Exercise: Explore population structure in a rice SNP panel

Use `data/rice_diversity_panel_synthetic.vcf` (10 rice accessions, 25 SNPs).

## Tasks
1. Parse the VCF and print a table of SNP ID, position, nearby gene (from the `INFO` field), and overall
   ALT allele frequency across all 10 accessions.
2. Compute simple **observed heterozygosity** per SNP: fraction of accessions with genotype `0/1`. Which
   SNP has the highest heterozygosity?
3. Using the group-split method from the code-along, identify every SNP with a strong frequency
   difference (>0.5) between the first 5 and last 5 accessions. How many are there, and are they
   clustered near particular genes in the `INFO` field?
4. **Discussion:** in real rice genetics, the species has two major subspecies, *indica* and *japonica*,
   which are genetically quite distinct. If your 10 accessions actually represented 5 indica + 5 japonica
   varieties, would you expect MANY SNPs across the genome to show this kind of group split, or just a
   few? What does that imply about using a handful of SNPs vs. thousands to infer population structure
   reliably?
5. **Stretch:** compute each accession's total number of ALT alleles across all 25 SNPs (a crude
   "genetic distance from reference" score) and rank the 10 accessions from most to least reference-like.

## Answer key
```python
def parse_vcf(path):
    samples = None
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]; continue
            fields = line.strip().split("\t")
            genotypes = dict(zip(samples, fields[9:]))
            records.append({"id": fields[2], "pos": fields[1], "info": fields[7], "genotypes": genotypes})
    return samples, records

samples, records = parse_vcf("../../data/rice_diversity_panel_synthetic.vcf")

def heterozygosity(genotypes):
    het = sum(1 for gt in genotypes.values() if gt == "0/1")
    return het / len(genotypes)

het_table = sorted(((r["id"], heterozygosity(r["genotypes"])) for r in records), key=lambda x: -x[1])
print("Most heterozygous SNP:", het_table[0])

alt_counts = {s: 0 for s in samples}
for r in records:
    for s, gt in r["genotypes"].items():
        alt_counts[s] += gt.replace("|", "/").split("/").count("1")
ranked = sorted(alt_counts.items(), key=lambda x: x[1])
print("Most reference-like accession:", ranked[0])
print("Most divergent-from-reference accession:", ranked[-1])
```
Discussion: true, genome-wide subspecies structure (like indica vs. japonica) should show up as a
consistent signal across **many** SNPs genome-wide, not just one or two — a single differentiated SNP could
easily be noise or a locally selected allele, which is exactly why real population-structure analyses
(e.g., STRUCTURE, PCA on genome-wide SNPs) use thousands to millions of markers, not a handful.
