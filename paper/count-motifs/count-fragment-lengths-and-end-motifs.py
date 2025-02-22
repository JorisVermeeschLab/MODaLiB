from typing import Optional, List, Dict

import os
import sys
import tqdm
import random
import time
import json

import numpy as np
import pysam

from modalib.core import Sample, ReferenceNucMap, chr_id_to_chr_name


ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, '..', '..', 'step1-bam-processing', 'out_')


def count_end_motifs(filepath) -> np.ndarray:
    sample = Sample(filepath)
    freqs = np.zeros((256, 500), dtype=np.uint64)
    for frag in sample:
        if frag['mapq'] < 60:
            continue
        i = frag['end_motif_5p_int']
        j = frag['insert_size']
        if j < 500:
            freqs[i, j] += 1
    return freqs


gc_code = 'GC114127'
filepath = os.path.join(DATA_DIR, f'{gc_code}.modalib')
np.save(f'{gc_code}-endmotifs-lengths.npy', count_end_motifs(filepath))
