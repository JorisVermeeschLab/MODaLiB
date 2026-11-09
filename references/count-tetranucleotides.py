import sys
sys.path.append('..')
from modalibpy.end_motifs import *

from Bio import SeqIO


fasta_file = 'hg38.fa'

freqs = np.zeros(256, dtype=int)

ALLOWED_CHROMOSOMES = {f"chr{c}" for c in range(1, 23)}

with open(fasta_file, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        print(record.id)
        if record.id in ALLOWED_CHROMOSOMES:
            sequence = str(record.seq)
            for i in range(len(sequence) - 4):
                motif = sequence[i:i+4].upper()
                if 'N' in motif:
                    continue
                k = str_motif_to_uint8(motif)
                freqs[k] += 1
#freqs /= np.sum(freqs)

print(list(freqs))
