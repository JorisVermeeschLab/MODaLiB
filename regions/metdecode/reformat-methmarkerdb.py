import os

import numpy as np
import pandas as pd


filepath = 'Colorectal_cancer.SMART.dmr.merge.DMRanno.txt'

markers = {}
with open(filepath, 'r') as f:
    lines = f.readlines()
    header = lines[0]
    for line in lines[1:]:
        elements = line.rstrip().split('\t')
        chrom = elements[0]
        start = int(elements[1])
        end = int(elements[2])
        length = int(elements[3])
        gene_name = elements[11] + '_' + elements[9]
        p_value = float(elements[8])

        if (p_value > 1e-50) or (length < 5000):
            continue

        if chrom.startswith('chr'):
            chrom = chrom[3:]
        if chrom == 'X':
            chr_id = 23
        elif chrom == 'Y':
            chr_id = 24
        elif chrom == 'MT':
            chr_id = 25
        elif chrom.isnumeric():
            chr_id = int(chrom)
        else:
            raise ValueError(chrom)

        coord = 300000000 * chr_id + start

        markers[coord] = (chr_id, start, end, gene_name)

coords = np.asarray(list(markers.keys()), dtype=np.uint64)
coords = np.sort(coords)

with open('dmrs-methmarkerdb-colorectal.csv', 'w') as f:
    for i, coord in enumerate(coords):
        chr_id, start, end, gene_name = markers[int(coord)]
        f.write(f'{chr_id},{start},{end},{gene_name}\n')
