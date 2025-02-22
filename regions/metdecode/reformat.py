import os

import numpy as np
import pandas as pd


filepath = 'kidney_atlas_tims.tims.txt'

markers = {}
with open(filepath, 'r') as f:
    lines = f.readlines()
    header = lines[0]
    for line in lines[1:]:
        elements = line.rstrip().split('\t')
        chrom = elements[0]
        start = int(elements[1])
        end = int(elements[2])

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

        markers[coord] = (chr_id, start, end, line)

coords = np.asarray(list(markers.keys()), dtype=np.uint64)
coords = np.sort(coords)

with open('kidney_atlas_tims.sorted.tsv', 'w') as f:
    f.write(header)
    for coord in coords:
        _, _, _, line = markers[int(coord)]
        f.write(line)

with open('dmrs-kidney.csv', 'w') as f:
    for i, coord in enumerate(coords):
        chr_id, start, end, _ = markers[int(coord)]
        f.write(f'{chr_id},{start},{end},M{i + 1}\n')
