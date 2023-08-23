# -*- coding: utf-8 -*-
#
#  nucleosomes.py
#
#  Copyright 2023 Antoine Passemiers <antoine.passemiers@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

import os
from typing import List

import numpy as np
import scipy.special

from modalib.core import Sample


# Each key is a chromosome ID, and each value is a tuple
# of the form (index, number of 50kb bins in chromosome).
CHR_LENGTHS = {
    1: (0, 4945),
    2: (1, 4860),
    3: (2, 3991),
    4: (3, 3826),
    5: (4, 3618),
    6: (5, 3418),
    7: (6, 3177),
    8: (7, 2926),
    9: (8, 2806),
    10: (9, 2708),
    11: (10, 2690),
    12: (11, 2647),
    13: (12, 2283),
    14: (13, 2128),
    15: (14, 2007),
    16: (15, 1777),
    17: (16, 1576),
    18: (17, 1523),
    19: (18, 1277),
    20: (19, 1249),
    21: (20, 939),
    22: (21, 994),
    23: (22, 3099),
    24: (23, 1156),
    25: (24, 1),
    26: (25, 4)
}


def nucleosome_scores(arr: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    # left  : |  20 | 146 |  20 | ... |  20 |     |     |
    # mid   : |     | 146 |  20 | ... |  20 | 146 |     |
    # right : |     |     |  20 | ... |  20 | 146 |  20 |
    cs = np.cumsum(arr, dtype=float)
    left = cs[20:-166] - cs[:-186]
    mid = cs[166:-20] - cs[20:-166]
    right = cs[186:] - cs[166:-20]

    shift = 4
    cs[shift:-186+shift] = left / 40.
    cs[20+shift:-166+shift] -= mid / 146.
    cs[166+shift:-20+shift] += right / 40.
    scores = scipy.special.expit(cs * gamma)
    scores[:186] = 0
    scores[-186:] = 0
    return scores


def save_piled_up_fragment_counts(filepaths: List[str], output_filepath: str):
    ref_map = {chr_id: np.zeros(CHR_LENGTHS[chr_id][1] * 50000, dtype=np.uint16) for chr_id in CHR_LENGTHS.keys()}
    for filepath in filepaths:
        sample = Sample(filepath)
        for read in sample:
            if read.chr_id in ref_map:
                chromosome = ref_map[read.chr_id]
                j = read.start
                length = read.size
                if j + length < len(chromosome):
                    chromosome[j:j+length] += 1

    ref_map = {str(key): value for key, value in ref_map.items()}
    np.savez(output_filepath, **ref_map)


def create_reference_nucleosome_map(input_filepath: str, output_folder: str):
    ref_map = dict(np.load(input_filepath))

    for chr_name in ref_map.keys():

        print(f'Processing chr{chr_name}...')

        for gamma in [0.05, 0.1, 0.2, 0.3, 0.4]:
            ref_map[chr_name] = nucleosome_scores(ref_map[chr_name], gamma=gamma)
        ref_map[chr_name] = (ref_map[chr_name] >= 0.5).astype(int)
        # ref_map[chr_name][0] = 0
        ref_map[chr_name] = np.diff(ref_map[chr_name])
        starts = np.where(ref_map[chr_name] > 0)[0]
        ends = np.where(ref_map[chr_name] < 0)[0][1:]
        m = min(len(starts), len(ends))
        mids = 0.5 * (starts[:m] + ends[:m])
        assert np.all(starts[:m] <= ends[:m])
        ref_map[chr_name] = mids

        ref_map[chr_name] = ref_map[chr_name]
        np.save(os.path.join(output_folder, f'{chr_name}.npy'), ref_map[chr_name])
