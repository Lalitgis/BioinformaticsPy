# Session 5 — Exercise: GC content & codon usage across species

Use `data/DREB2A_plant_orthologs.fasta` (5 species).

## Tasks
1. For each of the 5 orthologs, compute GC content and print a sorted table (highest GC first).
2. Compute codon usage for the Arabidopsis and maize sequences separately. Find the **top 3** codons
   used in each. Are they the same codons?
3. Compute the molecular weight of each ortholog's DNA sequence.
4. **Discussion:** if you were about to codon-optimize a bacterial *Bt* insecticidal gene for expression in
   maize, why would you specifically want maize's (not Arabidopsis's) codon usage table?
5. **Stretch:** write a function `gc_at_positions(seq_str)` that computes GC content separately at the 1st,
   2nd, and 3rd codon positions (this is a real technique — 3rd-position GC is usually most variable due
   to "wobble" degeneracy). Compare 1st/2nd/3rd position GC% for the FLC gene.

## Answer key
```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight

records = list(SeqIO.parse("../../data/DREB2A_plant_orthologs.fasta", "fasta"))
gc_table = sorted(((r.id, gc_fraction(r.seq)*100) for r in records), key=lambda x: -x[1])
for name, gc in gc_table:
    print(f"{name:40s} {gc:5.1f}%")

def codon_usage(seq_str):
    usage = {}
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        usage[codon] = usage.get(codon, 0) + 1
    return usage

for r in records:
    if "Arabidopsis" in r.id or "Zea" in r.id:
        u = codon_usage(str(r.seq))
        top3 = sorted(u.items(), key=lambda kv: -kv[1])[:3]
        print(r.id, "top codons:", top3)

for r in records:
    print(r.id, "MW:", round(molecular_weight(r.seq)))

def gc_at_positions(seq_str):
    pos = {1: "", 2: "", 3: ""}
    for i in range(0, len(seq_str) - 2, 3):
        pos[1] += seq_str[i]
        pos[2] += seq_str[i+1]
        pos[3] += seq_str[i+2]
    return {p: (s.count("G")+s.count("C"))/len(s)*100 for p, s in pos.items()}

flc = SeqIO.read("../../data/AT5G10140_FLC_synthetic_cds.fasta", "fasta").seq
print("GC% by codon position (FLC):", gc_at_positions(str(flc)))
```
