# Session 11 — Exercise: Build and compare two tree methods

Use `data/Solanaceae_matK_barcodes.fasta` (tomato, potato, pepper, petunia, tobacco).

## Tasks
1. Run the code-along to align the 5 sequences and build both a UPGMA and an NJ tree.
2. Print the two Newick strings (`tree.format("newick")` or just read the `.nwk` files) and compare them
   by eye. Do tomato and potato group together in **both** trees?
3. Using `dist_matrix`, print the single **smallest** pairwise distance and the single **largest**. Which
   two species are most similar? Which two are most different? Does this match real Solanaceae taxonomy
   (tomato/potato = genus *Solanum*; petunia and tobacco are more distantly related genera)?
4. **Discussion:** UPGMA assumes every lineage evolves at the same rate ("molecular clock"). Can you think
   of a biological reason two plant lineages might actually evolve at different rates (e.g., generation
   time, mutation repair efficiency, exposure to mutagens/UV)?
5. **Stretch:** re-run the whole pipeline using the `DREB2A_plant_orthologs.fasta` dataset instead. Does the
   resulting tree group rice and maize (both grasses) together, separate from the dicots?

## Answer key
```python
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

alignment = AlignIO.read("matk_aligned.fasta", "fasta")
calc = DistanceCalculator("identity")
dm = calc.get_distance(alignment)

pairs = []
names = dm.names
for i in range(len(names)):
    for j in range(i):
        pairs.append((names[i], names[j], dm[i, j]))

pairs.sort(key=lambda x: x[2])
print("Most similar pair:", pairs[0])
print("Most different pair:", pairs[-1])
```
Expect tomato/potato (both *Solanum*) to show the smallest distance, and petunia or tobacco (different
genera, longer independent evolutionary history within Solanaceae) to be among the most distant.
