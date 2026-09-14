# Session 6 — Exercise: Design your own Entrez fetch pipeline

**Important:** if your classroom has internet, do this exercise live against real NCBI. If not, write and
review the code on paper / walk through it together — the goal is understanding the *pattern*, which you can
run for real the next time you have a connection.

## Tasks
1. Write a function `search_count(term)` that returns just the total number of matching NCBI records for a
   search term (use the `"Count"` field from `esearch`'s result).
2. Compare hit counts for `"Oryza sativa[Organism]"` vs `"Arabidopsis thaliana[Organism]"` in the
   `nucleotide` database — which has more sequences deposited, and why might that be (hint: rice is a major
   global crop with huge breeding-program sequencing effort)?
3. Extend `fetch_first_hit()` from the code-along so it fetches the top **3** hits (not just 1) for a term,
   and saves them all to one multi-record FASTA file using `SeqIO.write(records, filename, "fasta")`.
4. Pick a crop and a trait you're curious about (e.g., "Solanum lycopersicum fruit ripening", "Triticum
   aestivum drought"), write the search term, and predict how many hits you'll get before running it.

## Answer key
```python
from Bio import Entrez, SeqIO
Entrez.email = "your_email@yourschool.edu"

def search_count(term, db="nucleotide"):
    handle = Entrez.esearch(db=db, term=term, retmax=0)
    result = Entrez.read(handle)
    handle.close()
    return int(result["Count"])

print("Rice hits:", search_count("Oryza sativa[Organism]"))
print("Arabidopsis hits:", search_count("Arabidopsis thaliana[Organism]"))

def fetch_top_n(term, n=3, db="nucleotide"):
    handle = Entrez.esearch(db=db, term=term, retmax=n)
    ids = Entrez.read(handle)["IdList"]
    handle.close()
    records = []
    for uid in ids:
        h = Entrez.efetch(db=db, id=uid, rettype="gb", retmode="text")
        records.append(SeqIO.read(h, "genbank"))
        h.close()
    return records

recs = fetch_top_n("Zea mays DREB2A", n=3)
SeqIO.write(recs, "maize_dreb2a_top3.fasta", "fasta")
```
