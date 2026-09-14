# Session 20 — Capstone Exercise (the whole course, in one project)

## The assignment
Working individually or in pairs, produce **one script** and **one half-page written report** answering
a research question of your choice, using only tools and datasets from this course. You may reuse the
DREB2A question from the code-along, or choose your own, for example:

- "Do the 5 Solanaceae species' matK sequences form a tree consistent with known genus-level taxonomy,
  and how much sequence divergence separates tomato/potato from the rest?" (uses Sessions 3, 9, 11-12)
- "Do all 4 candidate NBS-LRR resistance proteins carry both a P-loop and multiple LRR repeats, and how
  does their isoelectric point vary?" (uses Sessions 3, 14)
- "Which genes in the drought RNA-seq panel are both strongly differentially expressed AND annotated to
  a drought-relevant GO term?" (uses Sessions 17-18)
- Propose your own — any question answerable with 2+ combined techniques from the course qualifies.

## Requirements for the report
1. **Question** (1 sentence)
2. **Methods** (2-3 sentences: what data, what Biopython/pandas tools)
3. **Results** (the actual numbers you computed — no hand-waving)
4. **Conclusion** (what it means biologically) **+ one honest limitation** of your analysis

## Presentation (last 20 minutes of the session)
Each student/pair presents their one-page report in ~3 minutes: question, key number, one figure if
you made one, and your conclusion.

## Grading rubric (suggested, for the tutor)
| Criterion | Points |
|---|---|
| Question is answerable with course tools and clearly stated | 20 |
| Code runs without errors and produces the claimed numbers | 30 |
| Results are correctly interpreted (numbers match the conclusion) | 30 |
| Limitation is honestly and specifically identified (not generic) | 20 |

## Where to go from here (share with students)
- Re-run Session 6's Entrez code for real once you have normal internet — rebuild every dataset in this
  course with live NCBI records for your own species of interest.
- Try a real crop genome browser (Ensembl Plants, Phytozome) to explore genes beyond this course's panel.
- Learn DESeq2 (R) or PyDESeq2 (Python) for statistically rigorous differential expression.
- Explore GWAS properly with real tools like GAPIT, TASSEL, or plink, and real germplasm panel data
  (e.g., the Rice Diversity Panel, or SoyBase for soybean).
- Congratulations — you now have a working mental model of the entire plant bioinformatics analysis
  pipeline, and a full course you can keep re-teaching.
