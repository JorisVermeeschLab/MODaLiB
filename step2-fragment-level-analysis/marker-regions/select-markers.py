# -*- coding: utf-8 -*-
# select-markers.py
# author: Antoine Passemiers

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from recede.marker_selection import *


ROOT = os.path.dirname(os.path.abspath(__file__))


data = np.load('consolidated.atlas.npz')

D = data['coverage']
chromosomes = data['chrom']
coords = data['st']
tissue_names = data['tissue_names']
n_methylated = data['methylated']

# Remove some tissues
blacklist = {'WBC', 'breast', 'colon', 'placenta'}
mask = np.asarray([(tissue_name not in blacklist) for tissue_name in tissue_names], dtype=bool)
D = D[mask, :]
tissue_names = tissue_names[mask]
n_methylated = n_methylated[mask, :]

print(tissue_names)

# Ensure data is sorted by genomic coordinates
print('Sorting by genomic coordinates...')
values = chromosomes * 1e+10 + coords
idx = np.argsort(values)
chromosomes = chromosomes[idx]
coords = coords[idx]
D = D[:, idx]
n_methylated = n_methylated[:, idx]
for c in np.unique(chromosomes):
    mask = (chromosomes == c)
    _coords = coords[mask]
    diff = _coords[1:] - _coords[:-1]
    assert np.all(diff > 0)

D = np.maximum(D, 1)
M = n_methylated / D

print(tissue_names)
chrom, coords_start, coords_end, B, M = select_markers(M, D, chromosomes, coords,
        min_n_sites=25, max_distance=1000, tolerance=3, margin_threshold=0.4, n_iterations=3000000)

plt.rcParams['figure.figsize'] = [24, 12]
plt.subplot(1, 2, 1)
plt.imshow(M)
region_names = [f'chr{chrom[i]} {coords_start[i]}-{coords_end[i]}' for i in range(len(chrom))]
plt.xticks(range(len(region_names)), region_names, rotation=90)
plt.yticks(range(len(tissue_names)), tissue_names)
plt.subplot(1, 2, 2)
plt.imshow(B)
region_names = [f'chr{chrom[i]} {coords_start[i]}-{coords_end[i]}' for i in range(len(chrom))]
plt.xticks(range(len(region_names)), region_names, rotation=90)
plt.yticks(range(len(tissue_names)), tissue_names)
plt.tight_layout()
plt.savefig('atlas.png')

with open(os.path.join(ROOT, 'selected-markers.txt'), 'w') as f:
    for i in range(len(chrom)):
        f.write(f'chr{chrom[i]}, {coords_start[i]}, {coords_end[i]}\n')
