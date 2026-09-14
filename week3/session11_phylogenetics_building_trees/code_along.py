"""
Session 11 — Building phylogenetic trees (UPGMA & Neighbor-Joining)
"""
import subprocess
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

DATA = "../../data"

# Step 0: align first (same as Session 9)
subprocess.run([
    "clustalo", "-i", f"{DATA}/Solanaceae_matK_barcodes.fasta",
    "-o", "matk_aligned.fasta", "--outfmt=fasta", "--force",
], check=True)

alignment = AlignIO.read("matk_aligned.fasta", "fasta")
print(f"Aligned {len(alignment)} sequences, length {alignment.get_alignment_length()}")

# Step 1: distance matrix
calculator = DistanceCalculator("identity")
dist_matrix = calculator.get_distance(alignment)
print("\nDistance matrix:")
print(dist_matrix)

# Step 2: build trees two ways
constructor = DistanceTreeConstructor()

upgma_tree = constructor.upgma(dist_matrix)
print("\n=== UPGMA tree ===")
Phylo.draw_ascii(upgma_tree)

nj_tree = constructor.nj(dist_matrix)
print("\n=== Neighbor-Joining tree ===")
Phylo.draw_ascii(nj_tree)

# Save both for Session 12
Phylo.write(upgma_tree, "upgma_tree.nwk", "newick")
Phylo.write(nj_tree, "nj_tree.nwk", "newick")
print("\nSaved upgma_tree.nwk and nj_tree.nwk")
