import os
import time

import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from modalib.core import Sample, pairwise_kernels


# taskkill /IM python.exe /F
# python setup.py install
# python tmp.py

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, '..', 'step1-methyldackel', 'out')

# sample_names = ['GC121336', 'GC121345', 'GC114152']
sample_names = ['GC110680', 'GC118512']

samples = [Sample(f'../step1-methyldackel/out/{sample_name}') for sample_name in sample_names]

def load_hist(sample: Sample):
    gc_content = []
    for read in sample:
        if read.chr_id != 1:
            break
        gc_content.append(read.n_methylated / read.n_unmethylated)
    return gc_content


ys = load_hist(samples[0])
print(len(ys))
plt.hist(ys, density=True, alpha=0.4, bins=200, color='blue')
ys = load_hist(samples[1])
print(len(ys))
plt.hist(ys, density=True, alpha=0.4, bins=200, color='red')
plt.savefig('hist.png')
