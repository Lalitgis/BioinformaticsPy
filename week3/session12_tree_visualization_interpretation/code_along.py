"""
Session 12 — Tree visualization & interpretation
"""
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

DATA = "../../data"

subprocess.run(["clustalo", "-i", f"{DATA}/Solanaceae_matK_barcodes.fasta",
                 "-o", "matk_aligned.fasta", "--outfmt=fasta", "--force"], check=True)

alignment = AlignIO.read("matk_aligned.fasta", "fasta")
calc = DistanceCalculator("identity")
dm = calc.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Root using tobacco as outgroup
tree.root_with_outgroup({"name": "Nicotiana_tabacum_tobacco"})
tree.ladderize()   # sorts branches for a tidier figure

print("Rooted tree (ASCII):")
Phylo.draw_ascii(tree)

fig, ax = plt.subplots(figsize=(9, 4.5))
Phylo.draw(tree, axes=ax, do_show=False)
plt.tight_layout()
plt.savefig("solanaceae_tree.png", dpi=150)
print("\nSaved figure to solanaceae_tree.png")

# Print branch lengths for interpretation
print("\nTerminal branch lengths (evolutionary distance since divergence from its parent node):")
for clade in tree.get_terminals():
    print(f"  {clade.name:35s} {clade.branch_length:.4f}")
