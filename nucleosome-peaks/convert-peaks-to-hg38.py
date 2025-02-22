import tqdm
import liftover
import numpy as np


converter = liftover.get_lifter('hg19', 'hg38')

coords = []

with open('GSE71378_CA01.peaks.txt', 'r') as f:

    for line in tqdm.tqdm(f.readlines(), desc='Converting coordinates'):

        elements = line.rstrip().split('\t')
        if len(elements) == 2:
            chr_name, pos = elements

            # Convert from HG19 to HG38
            hits = converter.convert_coordinate(chr_name, int(pos))
            if len(hits) == 0:
                continue
            else:
                for chr_name, pos, _ in hits:
                    if '_' in chr_name:
                        continue
                    chr_id = chr_name.replace('chr', '')
                    if chr_id == 'X':
                        chr_id = 23
                    elif chr_id == 'Y':
                        chr_id = 24
                    else:
                        chr_id = int(chr_id)
                    coords.append([chr_id, pos])
coords = np.asarray(coords, dtype=np.ulonglong)

# Sorting coordinates
combined_coords = coords[:, 0] * 300_000_000 + coords[:, 1]
idx = np.argsort(combined_coords)
coords = coords[idx, :]


with open('GSE71378_CA01.peaks.hg38.txt', 'w') as f:
    for chr_id, pos in tqdm.tqdm(coords, desc='Writing'):
        chr_id = int(chr_id)
        if chr_id == 23:
            chr_name = 'chrX'
        elif chr_id == 24:
            chr_name = 'chrY'
        else:
            chr_name = f'chr{int(chr_id)}'
        pos = int(pos)
        f.write(f'{chr_name}\t{pos}\n')
