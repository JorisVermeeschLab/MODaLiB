import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from modalib.core import Sample

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(ROOT, '../step1-methyldackel/out')

"""
pos = np.load('nucleosomes/1.npy')

pos = pos[1:] - pos[:-1]
pos = pos[pos < 1000]

plt.hist(pos, bins=300)
plt.show()
"""

ys = []
filename = 'GC119350'
sample = Sample(os.path.join(DATA_PATH, filename))
for read in sample:
    ys.append(read.size)
    if len(ys) >= 100000:
        break
ys = np.asarray(ys)
ys = ys[ys < 1000]

plt.hist(ys, bins=300)
plt.savefig('ok.png')
