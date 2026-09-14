"""
Session 13 — Restriction sites & CAPS marker design
"""
from Bio import SeqIO
from Bio.Restriction import EcoRI, RestrictionBatch
from Bio.SeqUtils import MeltingTemp as mt

wt, mutant = list(SeqIO.parse("../../data/FLC_CAPS_marker_alleles.fasta", "fasta"))

print("--- EcoRI digestion simulation (CAPS marker) ---")
for allele in (wt, mutant):
    cuts = EcoRI.search(allele.seq)
    fragments = len(cuts) + 1
    print(f"{allele.id}")
    print(f"  EcoRI cut positions: {cuts}")
    print(f"  -> {fragments} fragment(s) expected on a gel\n")

print("--- Scanning with multiple enzymes at once ---")
rb = RestrictionBatch(["EcoRI", "BamHI", "HindIII", "PstI"])
result = rb.search(wt.seq)
for enzyme, positions in result.items():
    if positions:
        print(f"{enzyme}: cuts at {positions}")

print("\n--- Simple primer check ---")
def check_primer(primer_seq):
    tm = mt.Tm_Wallace(primer_seq)
    gc = (primer_seq.count("G") + primer_seq.count("C")) / len(primer_seq) * 100
    return {"length": len(primer_seq), "Tm": round(tm, 1), "GC%": round(gc, 1)}

forward_primer = str(wt.seq[:20])
reverse_primer = str(wt.seq[-20:].reverse_complement())

print("Forward primer:", forward_primer, check_primer(forward_primer))
print("Reverse primer:", reverse_primer, check_primer(reverse_primer))
