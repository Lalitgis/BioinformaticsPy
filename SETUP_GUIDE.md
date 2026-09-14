# Setup Guide — Plant Bioinformatics with Biopython

Do this once, before Session 1, on the machine(s) you'll teach with (your own laptop, or every machine in
a computer lab).

## 1. Install Python (3.9+)

Check first: `python3 --version`. If missing, install from [python.org](https://python.org) or via your
OS package manager.

## 2. Create an environment and install Python packages

**Option A — Anaconda/Miniconda (recommended for a teaching lab, easiest dependency management):**
```bash
conda create -n bioinfo python=3.11 -y
conda activate bioinfo
pip install biopython pandas matplotlib jupyter numpy
```

**Option B — plain venv + pip:**
```bash
python3 -m venv bioinfo-env
source bioinfo-env/bin/activate      # Windows: bioinfo-env\Scripts\activate
pip install biopython pandas matplotlib jupyter numpy
```

**Option C — Google Colab (zero local install, good if students lack laptops):**
```python
!pip install biopython
```
Datasets can be uploaded to Colab via the file-upload widget, or the whole `data/` folder can be mounted
from Google Drive.

## 3. Install command-line bioinformatics tools (needed from Week 2 onward)

These are standard, free, widely-used tools — not part of Biopython itself, but essential companions.

**Clustal Omega** (multiple sequence alignment — Sessions 9, 11, 12, 20):
```bash
# Ubuntu/Debian
sudo apt-get install clustalo
# macOS (Homebrew)
brew install clustal-omega
# Conda (any OS)
conda install -c bioconda clustalo
```

**NCBI BLAST+** (local sequence search — Session 10):
```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+
# macOS (Homebrew)
brew install blast
# Conda (any OS)
conda install -c bioconda blast
# Or download installers directly from NCBI:
# https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
```

## 4. Verify everything works

Run this check script (also usable as a first-day-of-class diagnostic for students):

```bash
python3 -c "
import Bio, pandas, matplotlib, numpy
print('Biopython:', Bio.__version__)
print('pandas:', pandas.__version__)
print('All Python packages OK')
"
clustalo --version
blastn -version
```

If all four commands print a version number with no errors, you're ready for Session 1.

## 5. About internet access and NCBI (Session 6)

Session 6 (Bio.Entrez) needs a normal internet connection to reach NCBI — this is different from every
other session, which works fully offline using the pre-built datasets in `data/`. If your classroom has
unreliable or filtered internet, you can still teach Session 6 conceptually from the slides and
`code_along.py` (which is correct, runnable code) without executing it live; students can try it later on
their own connection.

## 6. Folder layout reference

```
bioinfo-course/
├── README.html              <- open this first: course home, links to every session
├── SETUP_GUIDE.md            <- this file
├── data/                     <- every dataset, shared across sessions
├── week1/ ... week4/
│   └── sessionNN_slug/
│       ├── slides.html       <- open in a browser to present
│       ├── code_along.py     <- live-code this with your class
│       └── exercise.md       <- hands-on task + full answer key
```

## 7. A note on the datasets

Every file in `data/` is **synthetic-but-biologically-realistic**: built to mirror a real, named plant
gene, species, or published study, so the course works fully offline and reproducibly, with zero risk of a
broken internet connection derailing a 2-hour class. Each dataset's FASTA/VCF/CSV header documents which
real NCBI/TAIR record it represents and gives the exact search term to fetch the live version once you (or
your students) have normal internet access — Session 6 teaches exactly how to do that.

Good luck teaching — and enjoy watching students go from "what is a FASTA file" to running a full
comparative-genomics mini-pipeline in four weeks!
