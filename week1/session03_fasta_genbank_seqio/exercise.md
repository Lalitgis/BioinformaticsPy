# Session 3 — Exercise: Explore a plant gene family file

Use `data/DREB2A_plant_orthologs.fasta` (5 species) and `data/NBS_LRR_resistance_proteins.fasta`
(4 disease-resistance proteins) for this exercise.

## Tasks
1. Parse `DREB2A_plant_orthologs.fasta` and print a table: species id, sequence length, GC content
   (you may reuse your `gc_content()` function from Session 2, or Biopython's `Bio.SeqUtils.gc_fraction`).
2. Find which of the 5 orthologs is the **longest** and which is the **shortest**.
3. Parse `NBS_LRR_resistance_proteins.fasta` (a **protein** FASTA — same `SeqIO.parse` call works!) and
   print each protein's length.
4. Write a new FASTA file called `short_dreb2a.fasta` containing only the DREB2A orthologs under 460 bp.
5. **Discussion:** Why do orthologous genes (same gene, different species) have different lengths and
   sequences at all, if they do "the same job"? (Answer: molecular evolution — mutations accumulate at
   different rates in different lineages since their common ancestor; more distantly related species have
   had more time to diverge.)

## Answer key
```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

records = list(SeqIO.parse("../../data/DREB2A_plant_orthologs.fasta", "fasta"))
for r in records:
    print(f"{r.id:40s} len={len(r.seq):4d}  GC%={gc_fraction(r.seq)*100:.1f}")

lengths = {r.id: len(r.seq) for r in records}
longest = max(lengths, key=lengths.get)
shortest = min(lengths, key=lengths.get)
print("Longest:", longest, lengths[longest])
print("Shortest:", shortest, lengths[shortest])

proteins = list(SeqIO.parse("../../data/NBS_LRR_resistance_proteins.fasta", "fasta"))
for p in proteins:
    print(p.id, len(p.seq), "aa")

short_ones = [r for r in records if len(r.seq) < 460]
SeqIO.write(short_ones, "short_dreb2a.fasta", "fasta")
print(f"Wrote {len(short_ones)} records")
```
