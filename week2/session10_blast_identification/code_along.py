"""
Session 10 — Sequence identification with local BLAST
(remote NCBIWWW.qblast code is included below, commented, for classrooms with internet)
"""
import subprocess

DATA = "../../data"

# --- Step 1: build a local BLAST database from our plant gene family ---
subprocess.run([
    "makeblastdb", "-in", f"{DATA}/DREB2A_plant_orthologs.fasta",
    "-dbtype", "nucl", "-out", "dreb_db",
], check=True)

# --- Step 2: BLAST the mystery sequence against it ---
subprocess.run([
    "blastn",
    "-query", f"{DATA}/mystery_sequence_session10.fasta",
    "-db", "dreb_db",
    "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
    "-out", "blast_results.tsv",
], check=True)

print("BLAST results (query | subject | %identity | length | E-value | bitscore):")
with open("blast_results.tsv") as f:
    lines = f.readlines()
    for line in lines:
        print(" ", line.strip())

best_hit = lines[0].strip().split("\t")
print(f"\nBest hit: {best_hit[1]}  ({best_hit[2]}% identity, E-value={best_hit[4]})")
print("=> The mystery sequence is most likely a DREB2A-type ortholog closely related to this species.")

# --------------------------------------------------------------------------
# REMOTE BLAST AGAINST REAL NCBI (uncomment and run with internet access):
# --------------------------------------------------------------------------
# from Bio.Blast import NCBIWWW, NCBIXML
#
# with open(f"{DATA}/mystery_sequence_session10.fasta") as f:
#     query = f.read()
#
# result_handle = NCBIWWW.qblast("blastn", "nt", query)
# blast_record = NCBIXML.read(result_handle)
#
# for alignment in blast_record.alignments[:5]:
#     for hsp in alignment.hsps:
#         print(alignment.title)
#         print(f"  E-value: {hsp.expect}  Identity: {hsp.identities}/{hsp.align_length}")
