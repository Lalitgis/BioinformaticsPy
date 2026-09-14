# Session 2 — Exercise: Set up & explore your first gene

## Part A — Environment check (everyone)
1. Install Biopython using one of the 3 methods from the slides.
2. Run `code_along.py` and confirm you see "Biopython version: ..." printed with no errors.

## Part B — Coding practice
Using only what we covered today, write a script that:
1. Loads `data/AT5G10140_FLC_synthetic_cds.fasta` with `SeqIO.read()`.
2. Prints the sequence length.
3. Counts how many times each of A, C, G, T appears (use a dictionary, **not** Biopython shortcuts yet).
4. Computes the GC content **by hand** (don't use a built-in GC function — we build understanding first,
   shortcuts come in Session 5).
5. Prints the first 10 and last 10 bases of the sequence.

## Stretch goal (optional)
Write a function `count_codons(seq_str)` that splits a coding sequence into codons (groups of 3) and
returns a dictionary of codon → count. Test it on the FLC sequence. (We will build on this in Session 5.)

## Answer key
```python
from Bio import SeqIO

record = SeqIO.read("../../data/AT5G10140_FLC_synthetic_cds.fasta", "fasta")
seq = str(record.seq)

print("Length:", len(seq))

counts = {"A": 0, "C": 0, "G": 0, "T": 0}
for base in seq:
    counts[base] += 1
print("Base counts:", counts)

gc = (counts["G"] + counts["C"]) / len(seq) * 100
print(f"GC content: {gc:.2f}%")

print("First 10:", seq[:10])
print("Last 10:", seq[-10:])

def count_codons(seq_str):
    codons = {}
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        codons[codon] = codons.get(codon, 0) + 1
    return codons

print("Number of distinct codons used:", len(count_codons(seq)))
```
