"""
Session 5 — GC content, molecular weight, Tm, and codon usage
"""
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils import MeltingTemp as mt

flc = SeqIO.read("../../data/AT5G10140_FLC_synthetic_cds.fasta", "fasta").seq

print("--- FLC gene stats ---")
print(f"GC content: {gc_fraction(flc)*100:.2f}%")
print(f"Molecular weight: {molecular_weight(flc):.0f} g/mol")

short_primer_region = flc[:20]
print(f"Tm of first 20 bases (Wallace rule): {mt.Tm_Wallace(short_primer_region):.1f} C")

def codon_usage(seq_str):
    usage = {}
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        usage[codon] = usage.get(codon, 0) + 1
    return usage

usage = codon_usage(str(flc))
top5 = sorted(usage.items(), key=lambda kv: -kv[1])[:5]
print("Top 5 most-used codons in FLC:", top5)

print("\n--- Comparing GC content across all 5 DREB2A orthologs ---")
for rec in SeqIO.parse("../../data/DREB2A_plant_orthologs.fasta", "fasta"):
    gc = gc_fraction(rec.seq) * 100
    print(f"{rec.id:40s} GC={gc:5.1f}%")
