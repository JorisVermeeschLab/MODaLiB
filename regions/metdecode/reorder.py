import os

import numpy as np
import pandas as pd


filepath = 'v7_TCGA_regions_4CpG_250l_500bpdist.txt'
df = pd.read_csv(filepath, sep='\t', names=['cell_type', 'chrom', 'start', 'end', 'e1', 'e2'])

print(df)

markers = {}
for chrom, start, end, cell_type in zip(df['chrom'], df['start'], df['end'], df['cell_type']):
    start = int(start)
    end = int(end)
    cell_type = cell_type.replace('meth.', '')

    if cell_type in ['CESC', 'CEAD', 'OV']:
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

    markers[coord] = (chr_id, start, end, cell_type)

coords = np.asarray(list(markers.keys()), dtype=np.uint64)
coords = np.sort(coords)

with open('dmrs.csv', 'w') as f:
    for coord in coords:
        chr_id, start, end, cell_type = markers[int(coord)]
        f.write(f'{chr_id},{start},{end},{cell_type}\n')
