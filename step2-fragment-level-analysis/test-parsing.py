import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from modalib.core import Sample

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(ROOT, '../step1-bam-processing/out')

for filename in os.listdir(DATA_PATH):
    if filename.endswith('.modalib'):
        sample = Sample(os.path.join(DATA_PATH, filename))
        fragment = sample[len(sample) - 1]
        print(filename, fragment['start_chr_id'])
