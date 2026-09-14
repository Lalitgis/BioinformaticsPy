"""
Session 3 — FASTA & GenBank parsing with SeqIO
"""
from Bio import SeqIO

DATA = "../../data"

# --- Single record ---
flc = SeqIO.read(f"{DATA}/AT5G10140_FLC_synthetic_cds.fasta", "fasta")
print("Single record:", flc.id, len(flc.seq), "bp")

# --- Multi-record FASTA ---
print("\n--- DREB2A orthologs across 5 plant species ---")
orthologs = list(SeqIO.parse(f"{DATA}/DREB2A_plant_orthologs.fasta", "fasta"))
for rec in orthologs:
    print(f"{rec.id:40s} {len(rec.seq):4d} bp   {rec.description[-45:]}")

# --- GenBank record: sequence + annotation ---
print("\n--- GenBank record ---")
gb = SeqIO.read(f"{DATA}/FLC_teaching_record.gb", "genbank")
print("ID:", gb.id)
print("Organism:", gb.annotations.get("organism"))
print("Molecule type:", gb.annotations.get("molecule_type"))
print("Number of features:", len(gb.features))
for feat in gb.features:
    print(" -", feat.type, feat.location, dict(feat.qualifiers).get("gene", ""))

# --- Filter & write out ---
long_ones = [r for r in orthologs if len(r.seq) > 440]
SeqIO.write(long_ones, "filtered_output.fasta", "fasta")
print(f"\nWrote {len(long_ones)} of {len(orthologs)} DREB2A orthologs to filtered_output.fasta")
