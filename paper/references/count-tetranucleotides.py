import sys
sys.path.append('..')
from modalibpy.end_motifs import *

from Bio import SeqIO


fasta_file = 'genome.fa'

freqs = np.zeros(256, dtype=int)

with open(fasta_file, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        print(record.id)
        if record.id == 'MT':
            sequence = str(record.seq)
            for i in range(len(sequence) - 4):
                motif = sequence[i:i+4].upper()
                if 'N' in motif:
                    continue
                k = str_motif_to_uint8(motif)
                freqs[k] += 1
#freqs /= np.sum(freqs)

print(list(freqs))
