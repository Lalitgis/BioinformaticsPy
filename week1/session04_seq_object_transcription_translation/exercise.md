# Session 4 — Exercise: Transcribe & translate two DREB2A orthologs

Use `data/DREB2A_plant_orthologs.fasta`.

## Tasks
1. Load the file and pick the **Arabidopsis** and **rice** DREB2A-like sequences.
2. For each: transcribe to mRNA, translate to protein (use `to_stop=True`).
3. Print the protein lengths for both. Are they exactly the same length? If not, why might that be, given
   these are orthologs of the same transcription-factor family?
4. Reverse-complement the Arabidopsis sequence and try translating *that* directly (no `to_stop`). What do
   you get, and why does this demonstrate why strand orientation matters when you download a real gene
   from NCBI?
5. **Challenge:** write a function `find_start_codon(seq)` that scans a sequence and returns the index of
   the first `ATG` it finds in-frame with position 0 (i.e. index % 3 == 0). Test it on all 5 orthologs.

## Answer key
```python
from Bio import SeqIO

records = {r.id: r.seq for r in SeqIO.parse("../../data/DREB2A_plant_orthologs.fasta", "fasta")}
at = [v for k, v in records.items() if "Arabidopsis" in k][0]
os_ = [v for k, v in records.items() if "Oryza" in k][0]

for name, seq in [("Arabidopsis", at), ("Rice", os_)]:
    protein = seq.translate(to_stop=True)
    print(name, "protein length:", len(protein))

# reverse complement + translate (wrong strand) usually gives a very different,
# often much shorter, garbled protein full of early stop codons
print("Wrong-strand translation:", at.reverse_complement().translate())

def find_start_codon(seq):
    s = str(seq)
    for i in range(0, len(s) - 2, 3):
        if s[i:i+3] == "ATG":
            return i
    return -1

for name, seq in records.items():
    print(name, "first in-frame ATG at index:", find_start_codon(seq))
```
