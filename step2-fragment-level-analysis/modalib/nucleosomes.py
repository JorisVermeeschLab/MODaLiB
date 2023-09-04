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


def save_piled_up_fragment_counts(filepaths: List[str], output_filepath: str):
    """
    Piles up the fragments stored in the input files.

    For each base, counts the number of fragments spanning that locus.
    Intuitively, the higher the count of a specific locus, the higher the chances of the presence
    of the nucleosome, since DNA is more likely to be cleaved in the linker sequences.

    In the output file, each chromosome is represented as a separate NumPy array.
    For example, array for chromosome 6 is accessed using key "6" in the .npz file.

    Args:
        filepaths: List of modalib file paths.
        output_filepath: Where to store the .npz output file.
    """
    ref_map = {chr_id: np.zeros(CHR_LENGTHS[chr_id][1] * 50000, dtype=np.uint16) for chr_id in CHR_LENGTHS.keys()}
    for filepath in filepaths:
        sample = Sample(filepath)
        for frag in sample:  # Go through all the fragments in the sample
            if frag.start_chr_id in ref_map:
                chromosome = ref_map[frag.start_chr_id]
                j = frag.start_pos
                length = frag.size
                if j + length < len(chromosome):
                    # Each base position spanned by the fragment is incremented by 1
                    chromosome[j:j+length] += 1

    # Save fragments counts to a binary file
    ref_map = {str(key): value for key, value in ref_map.items()}
    np.savez(output_filepath, **ref_map)


def smooth_nucleosome_scores(arr: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    """
    Smooth nucleosome scores by applying a sliding window of size 186 bp (1 nucleosome + 2 linkers).

    If the fragment counts are high in a contiguous region of 146 bp, then the counts are expected
    to be lower in the 20 bp flanking regions. Nucleosome scores capture that idea, by subtracting
    the mean of scores from the flanking regions from the means of scores from the central 146 bp region.

    The scores are then scaled by `gamma` and being applied a sigmoid function in order to bound
    the scores in the [0, 1] range.

    Args:
        arr: Chromosome, where each value is either the number of fragments spanning the
            base position, or the likelihood of the presence of a nucleosome at that
            base position.
        gamma: Scaling factor

    Returns:
        Refined version of `arr`.
    """

    # Compute sums of scores in central 146 bp region and flanking 20 bp regions
    # left  : |  20 | 146 |  20 | ... |  20 |     |     |
    # mid   : |     | 146 |  20 | ... |  20 | 146 |     |
    # right : |     |     |  20 | ... |  20 | 146 |  20 |
    scores = np.cumsum(arr, dtype=float)
    left = scores[20:-166] - scores[:-186]
    mid = scores[166:-20] - scores[20:-166]
    right = scores[186:] - scores[166:-20]

    # Compute nucleosome score of each base position as the mean the central region
    # minus the means of the flanking regions
    shift = 4
    scores[shift:-186 + shift] = -left / 40.
    scores[20 + shift:-166 + shift] += mid / 146.
    scores[166 + shift:-20 + shift] -= right / 40.

    # Re-scale the scores and bound them in the [0, 1] range
    scores *= gamma
    scipy.special.expit(scores, out=scores)
    scores[:186] = 0
    scores[-186:] = 0
    return scores


def create_reference_nucleosome_map(input_filepath: str, output_folder: str, gamma0: float = 0.05):
    """
    Args:
        input_filepath: Path to the .npz containing the piled-up fragment counts.
        output_folder: Where to store the reference nucleosome map.
        gamma0: Scaling factor for the initial fragment counts.
            The lower, the lower the chances of getting a binary signal too early.
            The higher, the lower the chances of getting a smooth signal.
    """

    ref_map = dict(np.load(input_filepath))

    for chr_name in ref_map.keys():

        print(f'Processing chr{chr_name}...')

        # Iteratively smooth the nucleosome scores to obtain a binary signal
        for gamma in [gamma0, 0.1, 0.2, 0.3, 0.4]:
            ref_map[chr_name] = smooth_nucleosome_scores(ref_map[chr_name], gamma=gamma)

        # At this point, the signal should like binary, and alternate between regions of 1s
        # and regions of 0s. The center of each region of 1s is called as the center position of
        # a nucleosome.
        ref_map[chr_name] = (ref_map[chr_name] >= 0.5).astype(int)
        # ref_map[chr_name][0] = 0
        ref_map[chr_name] = np.diff(ref_map[chr_name])
        starts = np.where(ref_map[chr_name] > 0)[0]
        ends = np.where(ref_map[chr_name] < 0)[0][1:]
        m = min(len(starts), len(ends))
        mids = 0.5 * (starts[:m] + ends[:m])
        assert np.all(starts[:m] <= ends[:m])
        ref_map[chr_name] = mids

        # Save reference nucleosome map
        ref_map[chr_name] = ref_map[chr_name]
        np.save(os.path.join(output_folder, f'{chr_name}.npy'), ref_map[chr_name])
