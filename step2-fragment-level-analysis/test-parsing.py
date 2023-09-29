import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from modalib.core import Sample

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(ROOT, '../step1-bam-processing/out')


hist = {chr_id: 0 for chr_id in range(1, 27)}
chr_mismatches = 0
filename = 'GC110685.modalib'
sample = Sample(os.path.join(DATA_PATH, filename))
for frag in sample:
    if frag['start_chr_id'] != frag['end_chr_id']:
        chr_mismatches += 1
        print(f'Found chromosome mismatch: {frag["start_chr_id"]}-{frag["start_pos"]} != {frag["end_chr_id"]}-{frag["end_pos"]}')
    hist[frag['start_chr_id']] += 1

print(hist)
