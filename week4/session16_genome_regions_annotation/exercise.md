# Session 16 — Exercise: Build a mini genome annotation report

Use `data/FLC_teaching_record.gb`.

## Tasks
1. Run the code-along to generate `region_annotation.gff3` from the GenBank record.
2. Extend the GFF3 writer to also output the CDS feature (not just "gene") as its own GFF3 line.
3. Write a function `summarize_gff3(path)` that reads a GFF3 file and prints, for each feature type present
   (gene, CDS, ...), how many there are and their average length in bp.
4. **Discussion:** our teaching record has just 1 gene in ~700 bp (very dense!). Real Arabidopsis averages
   roughly 1 gene per ~4-5 kb; real maize (much larger, repeat-rich genome) averages roughly 1 gene per
   ~30-60 kb. If you were skimming a newly assembled 100 kb scaffold and found only 1 predicted gene, what
   two very different explanations might that have (hint: think "real biology" vs. "assembly/annotation
   problem")?
5. **Stretch:** simulate a slightly bigger toy region by taking 3 copies of the FLC record's gene feature,
   placing them 2,000 bp apart along an imaginary chromosome, and use `intergenic_distances()`-style logic
   to compute the gaps between them.

## Answer key
```python
def summarize_gff3(path):
    from collections import defaultdict
    counts = defaultdict(int)
    lengths = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            ftype = fields[2]
            length = int(fields[4]) - int(fields[3]) + 1
            counts[ftype] += 1
            lengths[ftype].append(length)
    for ftype in counts:
        avg = sum(lengths[ftype]) / len(lengths[ftype])
        print(f"{ftype}: {counts[ftype]} feature(s), average length {avg:.0f} bp")

summarize_gff3("region_annotation.gff3")
```
Discussion: a 100 kb scaffold with only 1 gene could mean (a) it genuinely sits in a real gene-poor,
repeat-rich region (common in large crop genomes), OR (b) the assembly/annotation pipeline failed to
predict genes there (fragmented assembly, low sequencing coverage, or a repeat region that confused the
gene-finder) — real genome projects follow up with additional evidence (RNA-seq mapping, homology to
related species) before concluding a region is truly gene-poor.
