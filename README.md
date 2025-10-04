# Project 1 — DNA Sequence Explorer

**Goal (kid-friendly):** Read little DNA strings, count letters, find start/stop, and translate to protein words.

## What you'll learn
- What a **FASTA** file looks like (`>name` line, then sequence lines)
- **GC%** (how many G and C letters)
- **ORFs** (start at `ATG`, stop at `TAA/TAG/TGA`)
- Translate DNA to amino acids using a **codon table**

## How to run
```bash
python src/analyze_sequences.py data/sequences.fasta
```
It prints a friendly summary and saves a plot to `results/gc_distribution.png`.

## Files
- `data/sequences.fasta` — tiny example DNA sequences
- `src/sequence_tools.py` — helper functions (GC%, ORFs, translation)
- `src/analyze_sequences.py` — command-line script
- `results/` — plots saved here
