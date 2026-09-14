# Session 7 — Exercise: Build a gene annotation report

Use `data/FLC_teaching_record.gb`.

## Tasks
1. Write a function `summarize_record(record)` (as sketched in the slides) that prints organism, total
   length, and every feature with its type, location, and length in bp.
2. Confirm computationally that `len(5'UTR) + len(CDS) + len(3'UTR) == len(whole record)`.
3. Extract the CDS feature, translate it, and print the **first and last 5 amino acids** of the protein.
4. **Discussion:** the CDS feature's `.qualifiers` dictionary already contains a pre-computed
   `"translation"`. Why would NCBI store the protein translation directly instead of making every user
   translate it themselves? (Think about: genetic code variants, RNA editing, computational conveniences
   for non-programmers browsing GenBank.)
5. **Stretch:** modify `summarize_record()` to also report the **percentage of the total record** each
   feature occupies (e.g., "CDS: 84.2% of record").

## Answer key
```python
from Bio import SeqIO

record = SeqIO.read("../../data/FLC_teaching_record.gb", "genbank")

def summarize_record(record):
    print(f"ID: {record.id}  Organism: {record.annotations.get('organism')}")
    print(f"Total length: {len(record.seq)} bp")
    total = len(record.seq)
    for feat in record.features:
        length = len(feat.location)
        pct = length / total * 100
        print(f"  {feat.type:10s} {feat.location}  {length} bp  ({pct:.1f}%)")

summarize_record(record)

lengths = {f.type: len(f.location) for f in record.features if f.type != "gene"}
print("Sum check:", sum(lengths.values()), "vs total", len(record.seq))

for feat in record.features:
    if feat.type == "CDS":
        protein = feat.extract(record.seq).translate(to_stop=True)
        print("First 5 aa:", protein[:5], "Last 5 aa:", protein[-5:])
```
