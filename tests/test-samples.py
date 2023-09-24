# -*- coding: utf-8 -*-
#
#  test-samples.py
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
import tempfile
from typing import Tuple, List

import numpy as np
from modalib.core import Fragment, Sample, pairwise_kernels, extract_features


ROOT = os.path.dirname(os.path.abspath(__file__))
NUC_FILEPATH = os.path.join(ROOT, 'nuc.ref.map.bed')


def create_sample(fragments: List[Fragment]) -> Sample:
    with tempfile.NamedTemporaryFile(delete=False) as f:
        for frag in fragments:
            byteorder = 'little'
            f.write(frag.start_chr_id.to_bytes(1, byteorder))
            f.write(frag.start_pos.to_bytes(4, byteorder))
            f.write(frag.end_chr_id.to_bytes(1, byteorder))
            f.write(frag.end_pos.to_bytes(4, byteorder))
            f.write(frag.size.to_bytes(2, byteorder))
            f.write(frag.n_methylated.to_bytes(2, byteorder))
            f.write(frag.n_unmethylated.to_bytes(2, byteorder))
            f.write(frag.snp_count.to_bytes(2, byteorder))
            f.write(frag.gc_count.to_bytes(2, byteorder))
            f.write(frag.s_count.to_bytes(2, byteorder))
            f.write(frag.n_inserted.to_bytes(2, byteorder))
            f.write(frag.n_deleted.to_bytes(2, byteorder))
        os.fsync(f)
    return Sample(f.name)


def test_iterator():
    fragments = [
        Fragment(start_chr_id=1, start_pos=17000, end_chr_id=1, end_pos=17501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=17000, end_chr_id=1, end_pos=17501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=17004, end_chr_id=1, end_pos=17501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=17005, end_chr_id=1, end_pos=17501, size=500, n_methylated=5, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=17006, end_chr_id=1, end_pos=17501, size=500, n_methylated=8, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=17007, end_chr_id=1, end_pos=17501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
    ]
    sample = create_sample(fragments)

    for i, frag in enumerate(sample):
        print(frag, fragments[i])
        assert frag == fragments[i]

    # Delete temporary files
    os.remove(sample.get_name())


def test_kernels():
    sample1 = create_sample([
        Fragment(start_chr_id=1, start_pos=16020, end_chr_id=1, end_pos=16521, size=500, n_methylated=4, n_unmethylated=16, snp_count=3, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=16020, end_chr_id=1, end_pos=16471, size=450, n_methylated=8, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=2, start_pos=18000, end_chr_id=2, end_pos=18071, size=70, n_methylated=0, n_unmethylated=16, snp_count=2, gc_count=36, s_count=140),
    ])
    sample2 = create_sample([
        Fragment(start_chr_id=1, start_pos=16000, end_chr_id=1, end_pos=16501, size=500, n_methylated=2, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=16030, end_chr_id=1, end_pos=16131, size=100, n_methylated=8, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=2, start_pos=18040, end_chr_id=2, end_pos=18111, size=70, n_methylated=12, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
    ])
    filepaths = [sample1.get_name(), sample2.get_name()]
    matrices = pairwise_kernels(filepaths, NUC_FILEPATH)
    assert np.any(matrices[1][0] != 0) and np.any(matrices[1][0] != 1)


test_kernels()
print('Finished')
