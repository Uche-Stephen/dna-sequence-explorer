# Minimal DNA utilities
CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',                
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
    'TGC':'C','TGT':'C','TGA':'*','TGG':'W',
}

STOP_CODONS = {'TAA','TAG','TGA'}

def read_fasta(path):
    seqs = {}
    name = None
    parts = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(parts).upper().replace(' ', '')
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line)
    if name:
        seqs[name] = ''.join(parts).upper().replace(' ', '')
    return seqs

def gc_content(seq:str)->float:
    seq = seq.upper()
    gc = sum(1 for c in seq if c in 'GC')
    atgc = sum(1 for c in seq if c in 'ATGC')
    return (gc/atgc*100.0) if atgc else 0.0

def translate(seq):
    aa = []
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)

def find_orfs(seq):
    # Return list of (start_index, end_index_exclusive, protein) in frame 0 only (simple).
    seq = seq.upper()
    orfs = []
    i = 0
    n = len(seq)
    while i < n-2:
        codon = seq[i:i+3]
        if codon == 'ATG':
            # walk until stop
            j = i
            while j < n-2:
                c = seq[j:j+3]
                if c in STOP_CODONS:
                    prot = translate(seq[i:j])
                    orfs.append((i, j+3, prot))
                    i = j + 3
                    break
                j += 3
            else:
                # no stop found
                i += 3
        else:
            i += 3
    return orfs
