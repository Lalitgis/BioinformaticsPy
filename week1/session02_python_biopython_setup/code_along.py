"""
Session 2 — Environment setup & first Biopython script
Run this after installing Biopython to confirm everything works.
"""
import Bio
from Bio.Seq import Seq
from Bio import SeqIO

print("Biopython version:", Bio.__version__)

# --- Python refresher: manual GC content on a plain string ---
def gc_content(seq_str):
    g = seq_str.count("G")
    c = seq_str.count("C")
    return round((g + c) / len(seq_str) * 100, 2)

demo = "ATGGCTAGCTAGGCGC"
print("Manual GC% of demo string:", gc_content(demo))

# --- Now the Biopython way ---
my_seq = Seq("ATGGCTAGCTAG")
print("\nSeq object:", my_seq)
print("Length:", len(my_seq))
print("Reverse complement:", my_seq.reverse_complement())

# --- Load today's real dataset ---
record = SeqIO.read("../../data/AT5G10140_FLC_synthetic_cds.fasta", "fasta")
print("\nLoaded record:", record.id)
print("Description:", record.description[:80], "...")
print("Sequence length:", len(record.seq))
print("First 30 bases:", record.seq[:30])
