# -*- coding: utf-8 -*-
# process-bed-files.py
# author: Antoine Passemiers

import os


ROOT = os.path.dirname(os.path.abspath(__file__))

IN_FOLDER = os.path.join(
    ROOT, '..', '..', '..', '..', 'methylation_data',
    'Methylation_atlas_replicates', 'filtered_bed')
OUT_FOLDER = os.path.join(ROOT, '..', 'step1-methyldackel', 'out', 'gdna')
if not os.path.isdir(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)

for filename in os.listdir(IN_FOLDER):
    with open(os.path.join(OUT_FOLDER, f'{filename}.out'), 'w') as fout:
        with open(os.path.join(IN_FOLDER, filename), 'r') as fin:
            print(f'Processing file {filename}...')
            for line in fin:
                elements = line.rstrip().split('\t')
                if len(elements) < 2:
                    continue
                chr_name = elements[0]
                start = int(elements[1])
                dmeth = int(elements[3])
                meth_ratio = float(elements[4])
                n_meth = int(round(meth_ratio * dmeth))
                n_unmeth = dmeth - n_meth
                fout.write(f'{chr_name}\t{start}\t{n_meth}\t{n_unmeth}\t0\t0\t2\t147\n')

print('Finished')
