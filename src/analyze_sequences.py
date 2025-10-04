import sys, os
import csv
import matplotlib.pyplot as plt
from sequence_tools import read_fasta, gc_content, find_orfs

# ===== marker to confirm you're running the right file =====
print("ANALYZE VERSION: v2")

# Make paths robust (works whether you run from root or from src)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')

def main(path):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Resolve input path relative to project root if it's not absolute
    if not os.path.isabs(path):
        candidate = os.path.join(PROJECT_ROOT, path)
        if os.path.exists(candidate):
            path = candidate

    seqs = read_fasta(path)
    print(f"Loaded {len(seqs)} sequences\n")

    gcs = []

    # Prepare summary.csv with header
    summary_path = os.path.join(RESULTS_DIR, 'summary.csv')
    with open(summary_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sequence_name', 'length_bp', 'gc_percent',
                         'orf_count', 'longest_orf_bp', 'longest_protein_aa'])

    for name, seq in seqs.items():
        gc = gc_content(seq)
        gcs.append(gc)
        orfs = find_orfs(seq)

        print(f"[{name}] length={len(seq)} bp | GC={gc:.1f}% | ORFs found={len(orfs)}")
        for k, (s, e, prot) in enumerate(orfs, start=1):
            print(f"  ORF{k}: {s}-{e} (len {e-s} bp) -> protein length {len(prot)} aa")

        # Longest ORF
        longest_orf_bp = 0
        longest_protein_aa = 0
        for (s, e, prot) in orfs:
            if (e - s) > longest_orf_bp:
                longest_orf_bp = (e - s)
                longest_protein_aa = len(prot)

        # Append row
        with open(summary_path, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([name, len(seq), f"{gc:.1f}", len(orfs),
                             longest_orf_bp, longest_protein_aa])
        print(f"  → wrote row to {summary_path}")

    # Plot GC%
    plt.figure()
    plt.hist(gcs, bins=10)
    plt.title('GC% distribution')
    plt.xlabel('GC%')
    plt.ylabel('Count')
    gc_plot = os.path.join(RESULTS_DIR, 'gc_distribution.png')
    plt.savefig(gc_plot, dpi=200, bbox_inches='tight')

    print(f"\nSaved plot -> {gc_plot}")
    print(f"Summary CSV -> {summary_path}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python src/analyze_sequences.py data/sequences.fasta")
        sys.exit(1)
    main(sys.argv[1])
