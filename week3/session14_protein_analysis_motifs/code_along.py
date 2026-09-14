"""
Session 14 — Protein analysis & motif finding
"""
import re
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

print("--- Physicochemical properties of each candidate R-protein ---")
records = list(SeqIO.parse("../../data/NBS_LRR_resistance_proteins.fasta", "fasta"))
for record in records:
    pa = ProteinAnalysis(str(record.seq))
    print(f"{record.id}")
    print(f"  Length: {len(record.seq)} aa")
    print(f"  MW: {pa.molecular_weight():.0f} Da")
    print(f"  Isoelectric point: {pa.isoelectric_point():.2f}")
    print(f"  Instability index: {pa.instability_index():.1f} "
          f"({'likely unstable' if pa.instability_index() > 40 else 'likely stable'})")

print("\n--- Scanning for the P-loop / Walker A motif (G-x(4)-GK[S/T]) ---")
pattern = re.compile(r"G.{4}GK[ST]")
for record in records:
    seq = str(record.seq)
    match = pattern.search(seq)
    if match:
        print(f"{record.id:30s} P-loop FOUND at position {match.start()}: {match.group()}")
    else:
        print(f"{record.id:30s} no P-loop motif found")
