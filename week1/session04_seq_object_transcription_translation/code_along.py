"""
Session 4 — Seq object: transcription & translation
"""
from Bio.Seq import Seq
from Bio import SeqIO

# --- Warm-up on a toy sequence ---
gene = Seq("ATGGCTAGCTAG")
print("Original       :", gene)
print("Complement     :", gene.complement())
print("Reverse comp.  :", gene.reverse_complement())
print("mRNA           :", gene.transcribe())
print("Protein (full) :", gene.translate())
print("Protein (clean):", gene.translate(to_stop=True))

# --- Now on the real FLC dataset ---
record = SeqIO.read("../../data/AT5G10140_FLC_synthetic_cds.fasta", "fasta")
seq = record.seq
print("\n--- FLC gene ---")
print("CDS length:", len(seq), "bp  (should be divisible by 3:", len(seq) % 3 == 0, ")")
print("Starts with ATG:", seq[:3] == "ATG")
print("Stop codon at end:", seq[-3:])

mrna = seq.transcribe()
print("mRNA (first 30):", mrna[:30])

protein = seq.translate(to_stop=True)
print("Protein length:", len(protein), "amino acids")
print("Protein sequence:\n", protein)

# --- Reading frame demo: what happens if we shift by 1 base? ---
print("\n--- Frame shift demo (why frame matters) ---")
print("Correct frame :", seq[:15].translate())
print("Shifted by +1 :", seq[1:16].translate())
print("Shifted by +2 :", seq[2:17].translate())
