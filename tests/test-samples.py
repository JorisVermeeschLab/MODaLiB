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

from modalib.core import Fragment, Sample


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
        os.fsync(f)
    return Sample(f.name)


def test_iterator():
    fragments = [
        Fragment(start_chr_id=1, start_pos=1000, end_chr_id=1, end_pos=1501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=1000, end_chr_id=1, end_pos=1501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=1004, end_chr_id=1, end_pos=1501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=1005, end_chr_id=1, end_pos=1501, size=500, n_methylated=5, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=1006, end_chr_id=1, end_pos=1501, size=500, n_methylated=8, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
        Fragment(start_chr_id=1, start_pos=1007, end_chr_id=1, end_pos=1501, size=500, n_methylated=4, n_unmethylated=16, snp_count=0, gc_count=36, s_count=140),
    ]
    sample = create_sample(fragments)

    for i, frag in enumerate(sample):
        assert frag == fragments[i]

    # Delete temporary files
    os.remove(sample.get_name())
