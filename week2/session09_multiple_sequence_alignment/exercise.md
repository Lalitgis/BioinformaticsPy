# Session 9 — Exercise: Align & mine the DREB2A gene family

Use `data/DREB2A_plant_orthologs.fasta` (5 species). Requires `clustalo` installed.

## Tasks
1. Run the alignment exactly as in the code-along and confirm all 5 sequences are now the same length.
2. Print, for each sequence, how many `-` gap characters were inserted (i.e., `str(record.seq).count("-")`).
   Which species needed the most gaps relative to the others?
3. Find the longest fully-conserved block (as in the code-along) and extract just that sub-alignment
   (all 5 sequences, only those columns). Translate it (careful: strip gaps first!) — what amino acids
   does this conserved core encode?
4. **Discussion:** if this conserved block turned out to sit inside the DNA-binding domain of the DREB2A
   transcription factor, what would that suggest about engineering drought tolerance by editing this gene
   — would you want to target this block with something like CRISPR, or avoid it? Why?
5. **Stretch:** compute a simple "conservation score" per column = (most common base's count) / (number of
   sequences), and plot it (a simple `print`-based bar using `"#" * int(score*20)` is fine if you don't have
   matplotlib handy) across the whole alignment.

## Answer key
```python
from Bio import AlignIO
import subprocess

subprocess.run(["clustalo", "-i", "../../data/DREB2A_plant_orthologs.fasta",
                 "-o", "aligned.fasta", "--outfmt=fasta", "--force"], check=True)
alignment = AlignIO.read("aligned.fasta", "fasta")

for record in alignment:
    print(record.id, "gaps:", str(record.seq).count("-"))

conserved_cols = [c for c in range(alignment.get_alignment_length())
                   if len(set(alignment[:, c])) == 1 and "-" not in alignment[:, c]]

best_run, cur_run, cur_start, best_start = 0, 0, 0, 0
for i, c in enumerate(conserved_cols):
    if i > 0 and c == conserved_cols[i-1] + 1:
        cur_run += 1
    else:
        cur_run, cur_start = 1, c
    if cur_run > best_run:
        best_run, best_start = cur_run, cur_start

block = alignment[:, best_start:best_start+best_run]
print("Conserved block sequence (species 1):", block[0].seq)

for col in range(0, alignment.get_alignment_length(), 10):
    column = alignment[:, col]
    from collections import Counter
    top_count = Counter(column).most_common(1)[0][1]
    score = top_count / len(alignment)
    print(f"col {col:4d}: {'#' * int(score*20)}")
```
