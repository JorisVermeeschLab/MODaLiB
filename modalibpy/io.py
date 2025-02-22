import os
from typing import Dict, List

import numpy as np


def load_profiles(filepaths: List[str]) -> Dict[str, np.ndarray]:
    if len(filepaths) == 0:
        return {}

    not_found = []
    for filepath in filepaths:
        if not os.path.exists(filepath):
            not_found.append(filepath)
    if not_found:
        raise FileNotFoundError(str(not_found))

    data = [np.load(filepath) for filepath in filepaths]
    keys = list(data[0].keys())
    out = {}
    for key in keys:

        invalid = []
        for filepath, d in zip(filepaths, data):
            if key not in d:
                invalid.append(filepath)
        #if invalid:
        #    raise KeyError(f'Key "{key}" missing in {invalid}')
        if not invalid:
            try:
                out[key] = [d[key] for d in data]
                out[key] = np.asarray(out[key])
            except ValueError:
                pass
    return out


def load_wig(filepath: str) -> Dict[int, np.ndarray]:
    values = {chr_id: [] for chr_id in range(1, 25)}
    chr_id = 1
    with open(filepath, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.rstrip()
        if len(line) < 2:
            continue
        if line.startswith('fixedStep'):
            chr_id = line.split(' ')[1].replace('chrom=chr', '')
            if chr_id == 'X':
                chr_id = 23
            elif chr_id == 'Y':
                chr_id = 24
            else:
                chr_id = float(chr_id)
        else:
            values[chr_id].append(float(line))
    for chr_id in range(1, 23):
        values[chr_id] = np.asarray(values[chr_id])
    return values
