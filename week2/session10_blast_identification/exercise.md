# Session 10 — Exercise: Identify three unknown plant sequences

## Setup
```bash
makeblastdb -in ../../data/DREB2A_plant_orthologs.fasta -dbtype nucl -out dreb_db
makeblastdb -in ../../data/Solanaceae_matK_barcodes.fasta -dbtype nucl -out matk_db
```

## Tasks
1. BLAST `data/mystery_sequence_session10.fasta` against `dreb_db` (as in the code-along). Report the best
   hit, its % identity, and its E-value in a one-sentence conclusion.
2. Pick any **one** sequence out of `data/Solanaceae_matK_barcodes.fasta` yourself, save it to its own
   FASTA file, then BLAST it against `matk_db` (which contains all 5, including itself). What % identity
   do you get for the "self-hit", and what does that tell you about a sanity check you should always run
   after building a new BLAST database?
3. Try BLASTing the DREB2A mystery sequence against `matk_db` instead (the wrong database on purpose).
   What happens to the E-values? Explain, in your own words, why a bad/no real match gives you weak
   E-values instead of an error.
4. **Discussion:** A colleague found a "significant" BLAST hit with 40% identity over only 15 bases. Should
   they trust it? What additional output field would settle the question?
5. **Stretch (needs internet):** uncomment the `NCBIWWW.qblast` block in the code-along, and remote-BLAST
   the mystery sequence against real NCBI's `nt` database. Compare the top real-world hit's description to
   what you found locally.

## Answer key
```python
import subprocess

subprocess.run(["makeblastdb", "-in", "../../data/Solanaceae_matK_barcodes.fasta",
                 "-dbtype", "nucl", "-out", "matk_db"], check=True)

# Self-hit sanity check
subprocess.run(["blastn", "-query", "../../data/Solanaceae_matK_barcodes.fasta",
                 "-db", "matk_db", "-outfmt", "6 qseqid sseqid pident evalue",
                 "-out", "self_check.tsv"], check=True)
with open("self_check.tsv") as f:
    for line in f:
        print(line.strip())
# Expect each sequence's top hit to be itself at 100% identity, E-value ~0 —
# this confirms the database was built correctly before trusting any other result.
```
Wrong-database BLAST (DREB2A query vs matK db) should return either no hits, or hits with poor E-values
(large, e.g. > 0.01) and low identity over a short alignment length — a strong signal of "no real homology,
just chance similarity," not a real ortholog relationship.
