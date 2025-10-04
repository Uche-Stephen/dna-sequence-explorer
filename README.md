# DNA Sequence Explorer

A simple bioinformatics project that reads DNA sequences from FASTA files, calculates GC content, detects open reading frames (ORFs), and produces both a histogram plot and a summary table.

---

##  Why this matters
GC% and ORF detection are fundamental in genome analysis. They are often the first steps in genome annotation, comparative genomics, and understanding coding potential in DNA sequences.

---

##  Features
- Parse FASTA files
- Calculate sequence length and GC content
- Identify open reading frames (ORFs) using `ATG` start and stop codons
- Export:
  - **Histogram plot** of GC% (`results/gc_distribution.png`)
  - **Summary table** with GC%, ORFs, and longest ORF (`results/summary.csv`)

---

## How to Run
Clone the repo and set up:

```bash
# Create environment
python -m venv .venv
.\.venv\Scripts\activate  # (Windows)
# or
source .venv/bin/activate  # (Mac/Linux)
```
## Install requirements
```bash
pip install -r requirements.txt
```

## Run the analyzer on a FASTA file:
```bash
python src/analyze_sequences.py data/sequences.fasta
```

Example Output

## results/gc_distribution.png
Histogram showing GC% across all sequences.

## results/summary.csv
sequence_name,length_bp,gc_percent,orf_count,longest_orf_bp,longest_protein_aa
seq1_example,39,56.4,1,24,7
seq2_example,31,32.3,0,0,0
seq3_example,32,34.4,1,9,2
my_seq,17,41.2,0,0,0

## What I learned
How to parse FASTA files programmatically
How to calculate GC% and detect ORFs
How to generate and save plots/CSV summaries in Python
Importance of documentation and reproducible analysis

## Next Steps / Future Improvements

Handle larger genomic datasets efficiently (e.g., multi-MB FASTA files)
Add motif search (e.g., find transcription factor binding sites)
Integrate with BLAST for sequence similarity search
Add error handling for malformed FASTA input
Explore codon usage statistics or GC-skew plots
Package as a command-line tool (pip install dna-explorer style)