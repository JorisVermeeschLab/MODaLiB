# -*- coding: utf-8 -*-
#
# core.pyx
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=True
#
#  Copyright 2022 Antoine Passemiers <antoine.passemiers@gmail.com>
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

from core cimport *

import os
import collections
from typing import Dict, List, Optional, Set

from libcpp.string cimport string
from cython.operator cimport dereference, preincrement
from cython.parallel import prange

import numpy as np
cimport numpy as cnp
cnp.import_array()


MIN_CHR_ID = 1
MAX_CHR_ID = 26



def coords_to_partition_(chr_id: int, pos: int, n_partitions: int) -> int:
    return coords_to_partition(chr_id, pos, n_partitions)


# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


Fragment = collections.namedtuple(
    'Fragment',
    ['chr_id', 'start', 'end', 'size', 'n_methylated', 'n_unmethylated', 'gc_count', 's_count']
)


cdef class Sample:

    cdef string filepath
    cdef _Sample sample_

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f'Could not open file {filepath}')
        self.filepath = <string>filepath.encode('utf-8')
        self.sample_ = _Sample(self.filepath)

    cdef _Sample get(self):
        return self.sample_

    def __iter__(self):
        while not self.sample_.eof():
            frag = self.sample_.get()
            yield Fragment(
                chr_id=frag.start_chr_id,
                start=frag.start_pos,
                end=frag.end_pos,
                size=frag.size,
                n_methylated=frag.n_methylated,
                n_unmethylated=frag.n_unmethylated,
                gc_count=frag.gc_count,
                s_count=frag.s_count
            )
            self.sample_.next()


cdef cnp.ndarray to_array(KernelMatrix &km):
    cdef size_t n = km.n_rows()
    cdef size_t m = km.n_cols()
    cdef size_t i, j
    cdef cnp.double_t[:, :] arr = np.empty((n, m), dtype=np.double)
    with nogil:
        for i in range(n):
            for j in range(m):
                arr[i, j] = km(i, j)
    return np.asarray(arr)


cdef vector[_Sample] samples_from_filepaths(filepaths: List[str]):
    cdef vector[_Sample] samples
    for filepath in filepaths:
        samples.push_back(_Sample(filepath))
    return samples


cdef unordered_map[cnp.uint8_t, vector[long]] load_ref_nuc_map(nuc_folder: str):
    # TODO: reference map is getting copied upon return -> optimize this
    cdef unordered_map[cnp.uint8_t, vector[long]] nuc_pos
    cdef long[:] pos
    cdef int i, chr_id
    for chr_id in range(1, 27):
        filepath = os.path.join(nuc_folder, f'{chr_id}.npy')
        if not os.path.exists(filepath):
            raise FileNotFoundError(
                f'Nucleosome position file (.npy format) is missing for chromosome "{chr_id}": '
                f'Could not locate file {filepath}.'
            )
        else:
            nuc_pos[chr_id] = vector[long]()
            pos = np.load(filepath).astype(int)
            with nogil:
                for i in range(pos.shape[0]):
                    nuc_pos[chr_id].push_back(pos[i])
    return nuc_pos


def pairwise_kernels(
        filepaths: List[str],
        nuc_folder: str
) -> Dict[int, List[np.ndarray]]:

    # Check number of samples
    n = len(filepaths)
    if n <= 1:
        raise ValueError(f'Number of samples should be at least 2. Got {n}.')

    # Load reference nucleosome map
    cdef unordered_map[cnp.uint8_t, vector[long]] nuc_pos = load_ref_nuc_map(nuc_folder)

    # Compute kernel matrices
    cdef unordered_map[cnp.uint8_t, vector[KernelMatrix]] matrices_ = compute_kernel_matrices(samples_from_filepaths(filepaths), nuc_pos)

    # Convert output kernel matrices to dict of lists of numpy arrays
    matrices = {}
    for chr_id in range(MIN_CHR_ID, MAX_CHR_ID + 1):
        if matrices_.find(chr_id) != matrices_.end():
            matrices[chr_id] = []
            for k in range(matrices_.at(chr_id).size()):
                matrices[chr_id].append(to_array(matrices_.at(chr_id).at(k)))

    return matrices


# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


cdef cnp.ndarray fm_to_array(FeatureMatrix &km):
    cdef size_t n = km.n_rows()
    cdef size_t m = km.n_cols()
    cdef size_t i, j
    cdef cnp.double_t[:, :] arr = np.empty((n, m), dtype=np.double)
    with nogil:
        for i in range(n):
            for j in range(m):
                arr[i, j] = km(i, j)
    return np.asarray(arr)


def extract_features(
        filepaths: List[str],
        nuc_folder: str
) -> Dict[int, List[np.ndarray]]:

    # Check number of samples
    n = len(filepaths)
    if n < 1:
        raise ValueError(f'Number of samples should be at least 1. Got {n}.')

    # Load reference nucleosome map
    cdef unordered_map[cnp.uint8_t, vector[long]] nuc_pos = load_ref_nuc_map(nuc_folder)

    # Compute feature matrices
    cdef unordered_map[cnp.uint8_t, vector[FeatureMatrix]] matrices_ = extract_all_features(samples_from_filepaths(filepaths), nuc_pos)

    # Convert output feature matrices to dict of list of numpy arrays
    matrices = {}
    for chr_id in range(MIN_CHR_ID, MAX_CHR_ID + 1):
        if matrices_.find(chr_id) != matrices_.end():
            matrices[chr_id] = []
            for k in range(matrices_.at(chr_id).size()):
                matrices[chr_id].append(fm_to_array(matrices_.at(chr_id).at(k)))

    return matrices
