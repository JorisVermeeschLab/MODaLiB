import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from modalib.core import Sample

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(ROOT, '../step1-bam-processing/out')

"""
pos = np.load('nucleosomes/1.npy')

pos = pos[1:] - pos[:-1]
pos = pos[pos < 1000]

plt.hist(pos, bins=300)
plt.show()
"""

xs, ys = [], []
filename = 'GC119350'
sample = Sample(os.path.join(DATA_PATH, filename))
for read in sample:
    ys.append(read.end - read.start)
    xs.append(read.size)
    if (len(ys) >= 100000):
        break
ys = np.asarray(ys)
xs = np.asarray(xs)
mask = np.logical_and(xs < 1000, ys < 1000)
ys = ys[mask]
xs = xs[mask]

plt.scatter(xs, ys, alpha=0.1)
plt.savefig('ok2.png')
