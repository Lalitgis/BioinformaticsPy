"""
Session 17 — Gene Ontology & functional annotation (manual, teaching version)
"""

# A small hand-curated GO-style annotation for our drought-study gene panel
go_annotation = {
    "DREB2A":   {"BP": "response to water deprivation", "MF": "DNA-binding transcription factor activity"},
    "RD29A":    {"BP": "response to water deprivation", "MF": "molecular function unknown"},
    "RD29B":    {"BP": "response to water deprivation", "MF": "molecular function unknown"},
    "NCED3":    {"BP": "abscisic acid biosynthetic process", "MF": "oxidoreductase activity"},
    "P5CS1":    {"BP": "proline biosynthetic process", "MF": "glutamate 5-kinase activity"},
    "RAB18":    {"BP": "response to abscisic acid", "MF": "molecular function unknown"},
    "RBCS1A":   {"BP": "photosynthesis, light reaction", "MF": "monooxygenase activity"},
    "RBCS2B":   {"BP": "photosynthesis, light reaction", "MF": "monooxygenase activity"},
    "CAB1":     {"BP": "photosynthesis, light harvesting", "MF": "chlorophyll binding"},
    "PSBA":     {"BP": "photosynthesis, light reaction", "MF": "electron transporter activity"},
    "ACTIN2":   {"BP": "cytoskeleton organization", "MF": "structural constituent of cytoskeleton"},
    "GAPDH":    {"BP": "glycolytic process", "MF": "oxidoreductase activity"},
    "UBQ10":    {"BP": "protein ubiquitination", "MF": "protein tag"},
    "TUB4":     {"BP": "microtubule-based process", "MF": "structural constituent of cytoskeleton"},
    "LEA14":    {"BP": "response to water deprivation", "MF": "molecular function unknown"},
    "ABI5":     {"BP": "response to abscisic acid", "MF": "DNA-binding transcription factor activity"},
    "PIP2_1":   {"BP": "water transport", "MF": "water channel activity"},
    "EXPA1":    {"BP": "cell wall organization", "MF": "structural molecule activity"},
    "CYP707A1": {"BP": "abscisic acid catabolic process", "MF": "oxidoreductase activity"},
    "COR15A":   {"BP": "response to cold", "MF": "molecular function unknown"},
}

all_genes = list(go_annotation.keys())
# "Drought-upregulated" gene list (matches the genes designed as "up" in Session 18's dataset)
drought_up = ["DREB2A", "RD29A", "RD29B", "NCED3", "P5CS1", "RAB18", "LEA14", "ABI5", "CYP707A1", "COR15A"]

def enrichment_check(gene_list, term, all_genes):
    hits_in_list = sum(1 for g in gene_list if go_annotation.get(g, {}).get("BP") == term)
    hits_in_background = sum(1 for g in all_genes if go_annotation.get(g, {}).get("BP") == term)
    frac_list = hits_in_list / len(gene_list)
    frac_background = hits_in_background / len(all_genes)
    return hits_in_list, len(gene_list), frac_list, hits_in_background, len(all_genes), frac_background

terms_to_check = ["response to water deprivation", "response to abscisic acid",
                   "photosynthesis, light reaction"]

print(f"Drought-upregulated gene list (n={len(drought_up)}): {drought_up}\n")
for term in terms_to_check:
    hits, n_list, frac_l, bg_hits, n_bg, frac_b = enrichment_check(drought_up, term, all_genes)
    fold = frac_l / frac_b if frac_b > 0 else float("inf")
    print(f"'{term}':")
    print(f"  In drought-up list: {hits}/{n_list} = {frac_l*100:.0f}%")
    print(f"  In whole panel (background): {bg_hits}/{n_bg} = {frac_b*100:.0f}%")
    print(f"  Enrichment fold: {fold:.1f}x\n")
