"""
Session 8 — Pairwise alignment with Bio.Align.PairwiseAligner
"""
from Bio import SeqIO, Align

records = {r.id: r.seq for r in SeqIO.parse("../../data/DREB2A_plant_orthologs.fasta", "fasta")}
at = [v for k, v in records.items() if "Arabidopsis" in k][0]
os_ = [v for k, v in records.items() if "Oryza" in k][0]

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

alignments = aligner.align(at, os_)
best = alignments[0]
print("Alignment score:", best.score)
print(best)

def percent_identity(alignment):
    a_str, b_str = str(alignment[0]), str(alignment[1])
    matches = sum(a == b for a, b in zip(a_str, b_str) if a != "-" and b != "-")
    return matches / len(a_str) * 100

print(f"Percent identity: {percent_identity(best):.1f}%")

print("\n--- Local alignment demo (same two sequences) ---")
aligner.mode = "local"
local_alignments = aligner.align(at, os_)
print("Best local alignment score:", local_alignments[0].score)
print(local_alignments[0])
