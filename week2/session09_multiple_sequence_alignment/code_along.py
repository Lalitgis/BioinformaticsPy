"""
Session 9 — Multiple sequence alignment with Clustal Omega + Bio.AlignIO
Requires clustalo installed and on PATH (see slides for install commands).
"""
import subprocess
from Bio import AlignIO

INPUT = "../../data/DREB2A_plant_orthologs.fasta"
OUTPUT = "dreb2a_aligned.fasta"

subprocess.run(
    ["clustalo", "-i", INPUT, "-o", OUTPUT, "--outfmt=fasta", "--force"],
    check=True,
)
print("Alignment written to", OUTPUT)

alignment = AlignIO.read(OUTPUT, "fasta")
print("\nNumber of sequences:", len(alignment))
print("Alignment length (with gaps):", alignment.get_alignment_length())

for record in alignment:
    print(f"{record.id:40s} {record.seq[:60]}...")

# --- Find fully conserved columns ---
conserved_cols = []
for col in range(alignment.get_alignment_length()):
    column = alignment[:, col]
    if len(set(column)) == 1 and "-" not in column:
        conserved_cols.append(col)

print(f"\nFully conserved columns: {len(conserved_cols)} out of {alignment.get_alignment_length()}")
print(f"Percent conserved: {len(conserved_cols)/alignment.get_alignment_length()*100:.1f}%")

# --- Find the longest run of consecutive conserved columns (a likely functional block) ---
best_run, cur_run, cur_start, best_start = 0, 0, 0, 0
for i, c in enumerate(conserved_cols):
    if i > 0 and c == conserved_cols[i-1] + 1:
        cur_run += 1
    else:
        cur_run = 1
        cur_start = c
    if cur_run > best_run:
        best_run = cur_run
        best_start = cur_start
print(f"Longest conserved block: {best_run} columns starting at alignment position {best_start}")
