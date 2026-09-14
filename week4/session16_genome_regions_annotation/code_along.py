"""
Session 16 — Genome regions, features & GFF3-style annotation
"""
from Bio import SeqIO

record = SeqIO.read("../../data/FLC_teaching_record.gb", "genbank")

# --- Treat our GenBank features as a mini "genome region" and summarize ---
genes = [f for f in record.features if f.type == "gene"]
print(f"Region length: {len(record.seq)} bp")
print(f"Number of genes in region: {len(genes)}")

total_gene_bp = sum(len(g.location) for g in genes)
density = total_gene_bp / len(record.seq) * 100
print(f"Gene density: {density:.1f}% of region is coding-gene territory")

# --- Write out a simple GFF3-style file from our GenBank record ---
with open("region_annotation.gff3", "w") as f:
    f.write("##gff-version 3\n")
    for feat in record.features:
        if feat.type == "gene":
            gene_id = feat.qualifiers.get("gene", ["unknown"])[0]
            start = int(feat.location.start) + 1  # GFF3 is 1-based
            end = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-"
            f.write(f"{record.id}\tteaching\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id}\n")

print("\nWrote region_annotation.gff3")
with open("region_annotation.gff3") as f:
    print(f.read())

# --- Simple GFF3 parser (for reading files like the one we just made) ---
def parse_gff3_line(line):
    fields = line.strip().split("\t")
    chrom, source, ftype, start, end, score, strand, phase, attrs = fields
    attr_dict = dict(pair.split("=") for pair in attrs.split(";") if "=" in pair)
    return {"chrom": chrom, "type": ftype, "start": int(start), "end": int(end),
            "strand": strand, "attrs": attr_dict}

with open("region_annotation.gff3") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parsed = parse_gff3_line(line)
        print(parsed)
