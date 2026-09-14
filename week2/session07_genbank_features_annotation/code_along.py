"""
Session 7 — GenBank feature tables & annotation
"""
from Bio import SeqIO

record = SeqIO.read("../../data/FLC_teaching_record.gb", "genbank")

print("ID:", record.id)
print("Organism:", record.annotations.get("organism"))
print("Total length:", len(record.seq), "bp")
print()

for feat in record.features:
    print(f"{feat.type:10s} {str(feat.location):20s} qualifiers={dict(feat.qualifiers)}")

print("\n--- Extracting and translating the CDS feature ---")
for feat in record.features:
    if feat.type == "CDS":
        cds_seq = feat.extract(record.seq)
        print("CDS length:", len(cds_seq))
        protein = cds_seq.translate(to_stop=True)
        print("Translated protein (first 40 aa):", protein[:40])
        # Compare to NCBI-provided translation stored in the qualifier, if present
        if "translation" in feat.qualifiers:
            ncbi_translation = feat.qualifiers["translation"][0]
            print("Matches stored translation qualifier:", str(protein) == ncbi_translation)

print("\n--- Extracting the 5'UTR and 3'UTR ---")
for feat in record.features:
    if feat.type in ("5\'UTR", "3\'UTR"):
        utr_seq = feat.extract(record.seq)
        print(feat.type, "length:", len(utr_seq), "sequence:", utr_seq)
