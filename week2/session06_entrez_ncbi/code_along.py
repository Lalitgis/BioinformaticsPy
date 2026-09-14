"""
Session 6 — Bio.Entrez (NCBI programmatic access)

NOTE: this sandbox has no outbound internet access to NCBI, so this script is
written to run correctly *as-is* on any normal machine with internet. Read it
top-to-bottom with students; if you have internet in class, run it live.
"""
from Bio import Entrez, SeqIO

Entrez.email = "your_email@yourschool.edu"   # <-- put a real email here

def fetch_first_hit(term, db="nucleotide"):
    """Search NCBI, return the first hit as a SeqRecord (GenBank format)."""
    search_handle = Entrez.esearch(db=db, term=term, retmax=1)
    search_result = Entrez.read(search_handle)
    search_handle.close()

    ids = search_result["IdList"]
    if not ids:
        print(f"No hits for: {term}")
        return None

    fetch_handle = Entrez.efetch(db=db, id=ids[0], rettype="gb", retmode="text")
    record = SeqIO.read(fetch_handle, "genbank")
    fetch_handle.close()
    return record

if __name__ == "__main__":
    # Example searches to try live once you have internet access:
    terms = [
        "Arabidopsis thaliana FLC mRNA complete cds",
        "Oryza sativa SD1 gene",
        "Zea mays DREB2A",
    ]
    for term in terms:
        print(f"\nSearching: {term}")
        rec = fetch_first_hit(term)
        if rec:
            print(" ->", rec.id, len(rec.seq), "bp,", rec.annotations.get("organism"))
            SeqIO.write(rec, f"{rec.id}.gb", "genbank")
            print("   saved to", f"{rec.id}.gb")
