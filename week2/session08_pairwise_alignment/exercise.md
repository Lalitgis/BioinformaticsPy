# Session 8 — Exercise: Align every pair of DREB2A orthologs

Use `data/DREB2A_plant_orthologs.fasta` (5 species).

## Tasks
1. Load all 5 sequences into a dictionary of `{name: seq}`.
2. Using `Align.PairwiseAligner` in **global** mode, compute the percent identity between every unique pair
   of the 5 sequences (10 pairs total). Print a sorted table from most-to-least similar.
3. Which pair is most similar? Which is least similar? Does this match what you'd predict from plant
   taxonomy (e.g., is maize more similar to rice — both grasses — than to Arabidopsis?)
4. Repeat the Arabidopsis-vs-rice comparison in **local** mode. Is the percent identity of the
   best local sub-alignment higher or lower than the global alignment? Why does that make sense?
5. **Discussion:** if two genes are only 60% identical at the DNA level, could their proteins still be
   nearly identical? (Hint: think about codon degeneracy from Session 4/5 — synonymous mutations don't
   change the protein.)

## Answer key
```python
from Bio import SeqIO, Align
from itertools import combinations

records = {r.id: r.seq for r in SeqIO.parse("../../data/DREB2A_plant_orthologs.fasta", "fasta")}

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

def percent_identity(alignment):
    a_str, b_str = str(alignment[0]), str(alignment[1])
    matches = sum(a == b for a, b in zip(a_str, b_str) if a != "-" and b != "-")
    return matches / len(a_str) * 100

results = []
for name1, name2 in combinations(records, 2):
    aln = aligner.align(records[name1], records[name2])[0]
    pid = percent_identity(aln)
    results.append((name1, name2, pid))

results.sort(key=lambda x: -x[2])
for n1, n2, pid in results:
    print(f"{n1:35s} vs {n2:35s} {pid:5.1f}%")
```
