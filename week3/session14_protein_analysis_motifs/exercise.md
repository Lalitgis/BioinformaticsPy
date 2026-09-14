# Session 14 — Exercise: Candidate resistance-gene screening

Use `data/NBS_LRR_resistance_proteins.fasta` (4 synthetic R-proteins).

## Tasks
1. Run the code-along. Confirm all 4 proteins carry the P-loop motif (they were designed to).
2. Write a second regex to search for the LRR-like repeat pattern `L..L.L..N.L` (leucine-rich repeats have
   a loose consensus like this). How many matches does each protein have? More matches suggests a longer
   LRR domain (more potential specificity for recognizing different pathogen effectors).
3. Compute and compare the isoelectric point (pI) of all 4 proteins. Proteins with pI below 7 are
   negatively charged at neutral pH; above 7, positively charged. Does this vary a lot across the 4?
4. **Discussion:** if you sequenced a new crop's genome and found 50 previously-unknown genes with a
   valid P-loop AND at least 3 LRR-like repeats, would you conclude they are definitely functional
   resistance genes? What would you still need to test experimentally to be sure?
5. **Stretch:** modify the motif regex to be *stricter* (e.g., require the exact classic P-loop
   `GMGGVGKT`) versus the looser `G.{4}GK[ST]`. Does the stricter pattern still match all 4 (it was built
   into all of them in the code-along)?

## Answer key
```python
import re
from Bio import SeqIO

lrr_pattern = re.compile(r"L..L.L..N.L")
records = list(SeqIO.parse("../../data/NBS_LRR_resistance_proteins.fasta", "fasta"))
for r in records:
    matches = lrr_pattern.findall(str(r.seq))
    print(r.id, "LRR-like repeat count:", len(matches))
```
Discussion answer: a motif match is only *suggestive* — true confirmation of resistance function requires
functional assays (e.g., transient expression + pathogen challenge, or a resistance/susceptibility
segregation test in a mapping population). Bioinformatics narrows thousands of genes down to a manageable
shortlist of candidates; it does not replace the wet-lab experiment.
