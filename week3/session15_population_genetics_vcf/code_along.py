"""
Session 15 — Population genetics & VCF parsing (pure Python, no extra library)
"""

def parse_vcf(path):
    samples = None
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue
            fields = line.strip().split("\t")
            chrom, pos, snp_id, ref, alt, qual, filt, info = fields[:8]
            genotypes = fields[9:]
            records.append({
                "chrom": chrom, "pos": int(pos), "id": snp_id,
                "ref": ref, "alt": alt, "info": info,
                "genotypes": dict(zip(samples, genotypes)),
            })
    return samples, records

def allele_frequency(genotypes):
    alt_count, total = 0, 0
    for gt in genotypes.values():
        alleles = gt.replace("|", "/").split("/")
        alt_count += alleles.count("1")
        total += len(alleles)
    return alt_count / total if total else 0.0

samples, records = parse_vcf("../../data/rice_diversity_panel_synthetic.vcf")
print(f"Loaded {len(records)} SNPs across {len(samples)} rice accessions")
print("Samples:", samples)

print("\n--- Allele frequency of first 8 SNPs ---")
for rec in records[:8]:
    freq = allele_frequency(rec["genotypes"])
    print(f"{rec['id']} (pos {rec['pos']}, near {rec['info']}): ALT freq = {freq:.2f}")

print("\n--- Checking for population structure (first 5 vs last 5 accessions) ---")
group_a = samples[:5]
group_b = samples[5:]
structured_snps = 0
for rec in records:
    fa = allele_frequency({s: rec["genotypes"][s] for s in group_a})
    fb = allele_frequency({s: rec["genotypes"][s] for s in group_b})
    if abs(fa - fb) > 0.5:
        structured_snps += 1
        print(f"{rec['id']}: group A={fa:.2f}  group B={fb:.2f}  (diff={abs(fa-fb):.2f})")

print(f"\n{structured_snps} of {len(records)} SNPs show strong group differentiation "
      f"(possible population structure signal)")
