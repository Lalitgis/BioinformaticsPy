# Session 13 — Exercise: Design a complete CAPS genotyping assay

Use `data/FLC_CAPS_marker_alleles.fasta` (wild-type and mutant FLC alleles).

## Tasks
1. Confirm computationally that the two alleles differ by exactly one base (hint: loop through both
   strings with `zip()` and count mismatches — they should be the same length).
2. Run the EcoRI digestion simulation from the code-along. Report clearly: how many bands would each
   allele show on a gel, and could you distinguish a **heterozygous** plant (one copy of each allele) if
   you ran a PCR product that was a 50/50 mixture of both alleles amplified together? (Hint: think about
   what bands would appear if you digested a mixture.)
3. Try 2 other enzymes from `RestrictionBatch` on both alleles. Do either of them also detect a
   difference between wild-type and mutant? (Most won't — CAPS markers work with a *specific* enzyme
   chosen precisely because it detects *that* SNP.)
4. Design a forward and reverse primer pair (20 bp each) that would amplify the entire allele sequence.
   Check both primers' Tm and GC% with `check_primer()`. Are their Tm values close enough to work well
   together in one PCR reaction (within ~3°C)?
5. **Discussion:** why is a CAPS marker often cheaper and faster for screening 500 breeding seedlings than
   sending all 500 samples for DNA sequencing?

## Answer key
```python
from Bio import SeqIO

wt, mutant = list(SeqIO.parse("../../data/FLC_CAPS_marker_alleles.fasta", "fasta"))
diffs = [i for i, (a, b) in enumerate(zip(str(wt.seq), str(mutant.seq))) if a != b]
print("Number of differing positions:", len(diffs), "at index/indices:", diffs)

# A heterozygous plant's PCR product is a mix of both allele copies. After digestion you'd see
# BOTH the cut fragments (from the wild-type copy) AND the full-length uncut fragment (from the
# mutant copy) on the same gel lane — a 3-band pattern distinguishing heterozygotes from either
# homozygote, which is exactly why CAPS markers are useful for co-dominant genotyping.
```
