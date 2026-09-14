# Session 12 — Exercise: Root, draw, and interpret

Use the tree built in Session 11 (or rebuild it) from `data/Solanaceae_matK_barcodes.fasta`.

## Tasks
1. Root the tree using tobacco as the outgroup, then draw both the ASCII tree and a saved PNG figure.
2. List all terminal branch lengths from longest to shortest. Which species has accumulated the most
   sequence change since it diverged from its nearest relative in this tree?
3. Try rooting with **pepper** instead of tobacco as the outgroup. Does the overall grouping of
   tomato+potato change? What *does* change?
4. Write a 2–3 sentence "briefing note" as if you were advising a tomato breeding program: which of the
   other 4 species in this tree would be the best candidate for introducing new genetic diversity via
   traditional crossing, based on relatedness alone (setting aside real crossing-compatibility biology)?
5. **Discussion:** the exercise says "setting aside real crossing compatibility." Why can't relatedness in
   a phylogenetic tree alone guarantee two species can actually be crossed to produce viable offspring?

## Answer key
```python
tree.root_with_outgroup({"name": "Capsicum_annuum_pepper"})
Phylo.draw_ascii(tree)
# Tomato and potato should still group together regardless of which outgroup you choose,
# because that grouping reflects genuine shared ancestry (both genus Solanum) — only the
# root position/outgroup placement changes, not the core topology among the remaining species.
```
Real biology note: reproductive isolation (chromosome number differences, hybrid sterility, pre-zygotic
barriers like pollen incompatibility) can block crossing between species that still look "close" on a
sequence-similarity tree — sequence similarity is necessary-ish evidence of relatedness, not sufficient
proof of cross-compatibility.
