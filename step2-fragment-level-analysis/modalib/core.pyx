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


Read = collections.namedtuple(
    'Read',
    ['chr_id', 'start', 'size', 'n_methylated', 'n_unmethylated', 'gc_count', 's_count']
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
            read = self.sample_.get()
            yield Read(
                chr_id=read.start_chr_id,
                start=read.start_pos,
                size=read.size,
                n_methylated=read.n_methylated,
                n_unmethylated=read.n_unmethylated,
                gc_count=read.gc_count,
                s_count=read.s_count
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


cdef vector[_Sample] samples_from_filepaths(vector[string]& filepaths) nogil:
    cdef vector[_Sample] samples
    cdef vector[string].iterator it = filepaths.begin()
    with nogil:
        while it != filepaths.end():
            samples.push_back(_Sample(dereference(it)))
            preincrement(it)
    return samples


cdef unordered_map[cnp.uint8_t, vector[KernelMatrix]] compute_kernel_matrices_whole_genome(
        vector[string]& filepaths,
        unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
) nogil:
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[bool] chr_mask
    chr_mask.push_back(False)
    for _ in range(26):
        chr_mask.push_back(True)
    return compute_kernel_matrices(samples, chr_mask, nuc_pos)


cdef unordered_map[cnp.uint8_t, vector[KernelMatrix]] compute_kernel_matrices_multiple_chromosomes(
        vector[string]& filepaths,
        set chr_ids,
        unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
):
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[bool] chr_mask
    for c in range(27):
        chr_mask.push_back(c in chr_ids)
    return compute_kernel_matrices(samples, chr_mask, nuc_pos)


cdef unordered_map[cnp.uint8_t, vector[KernelMatrix]] compute_kernel_matrices_one_chromosome(
        vector[string]& filepaths,
        cnp.uint8_t chr_id,
        unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
) nogil:
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[bool] chr_mask
    for c in range(27):
        chr_mask.push_back(c == chr_id)
    return compute_kernel_matrices(samples, chr_mask, nuc_pos)


def pairwise_kernels(
        samples: List[Sample],
        nuc_folder: str,
        chr_ids: Optional[Set[int]] = None,
        multithreaded: bool = False
) -> Dict[int, List[np.ndarray]]:

    n = len(samples)
    if n <= 1:
        raise ValueError(f'Number of samples should be at least 2. Got {n}.')

    cdef Sample sample
    cdef vector[_Sample] samples_
    cdef vector[string] filepaths
    cdef int i
    for i in range(len(samples)):
        sample = samples[i]
        samples_.push_back(sample.get())
        filepaths.push_back(sample.filepath)

    cdef unordered_map[cnp.uint8_t, vector[long]] nuc_pos
    cdef long[:] pos
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

    matrices = {}

    cdef KernelMatrix km
    cdef cnp.double_t[:, :] arr
    cdef unordered_map[cnp.uint8_t, vector[KernelMatrix]] matrices_, matrices_i_
    cdef int k, c
    if not multithreaded:
        if chr_ids is not None:
            matrices_ = compute_kernel_matrices_multiple_chromosomes(filepaths, chr_ids, nuc_pos)
        else:
            matrices_ = compute_kernel_matrices_whole_genome(filepaths, nuc_pos)
    else:
        # for c in range(1, 27):
        for c in prange(1, 27, nogil=True, num_threads=27, schedule='dynamic'):
            matrices_i_ = compute_kernel_matrices_one_chromosome(filepaths, c, nuc_pos)
            matrices_[c] = matrices_i_[c]
    for chr_id in range(1, 27):
        if matrices_.find(chr_id) != matrices_.end():
            matrices[chr_id] = []
            for k in range(matrices_[chr_id].size()):
                matrices[chr_id].append(to_array(matrices_[chr_id][k]))

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


cdef unordered_map[cnp.uint8_t, vector[FeatureMatrix]] compute_feature_matrices_whole_genome(
        vector[string]& filepaths,
        unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
) nogil:
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[bool] chr_mask
    chr_mask.push_back(False)
    for _ in range(27):
        chr_mask.push_back(True)
    return extract_all_features(samples, chr_mask, nuc_pos)


cdef unordered_map[cnp.uint8_t, vector[FeatureMatrix]] compute_feature_matrices_one_chromosome(
        vector[string]& filepaths,
        cnp.uint8_t chr_id,
        unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
) nogil:
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[bool] chr_mask
    for c in range(27):
        chr_mask.push_back(c == chr_id)
    return extract_all_features(samples, chr_mask, nuc_pos)


def extract_features(
        samples: List[Sample],
        nuc_folder: str,
        multithreaded: bool = False
) -> Dict[int, List[np.ndarray]]:

    n = len(samples)
    if n <= 1:
        raise ValueError(f'Number of samples should be at least 2. Got {n}.')

    cdef Sample sample
    cdef vector[_Sample] samples_
    cdef vector[string] filepaths
    cdef int i
    for i in range(len(samples)):
        sample = samples[i]
        samples_.push_back(sample.get())
        filepaths.push_back(sample.filepath)

    cdef unordered_map[cnp.uint8_t, vector[long]] nuc_pos
    cdef long[:] pos
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

    matrices = {}

    cdef FeatureMatrix km
    cdef cnp.double_t[:, :] arr
    cdef unordered_map[cnp.uint8_t, vector[FeatureMatrix]] matrices_, matrices_i_
    cdef int k, c
    if not multithreaded:
        matrices_ = compute_feature_matrices_whole_genome(filepaths, nuc_pos)
    else:
        # for c in range(1, 27):
        for c in prange(1, 27, nogil=True, num_threads=27, schedule='dynamic'):
            matrices_i_ = compute_feature_matrices_one_chromosome(filepaths, c, nuc_pos)
            matrices_[c] = matrices_i_[c]
    for chr_id in range(1, 27):
        if matrices_.find(chr_id) != matrices_.end():
            matrices[chr_id] = []
            for k in range(matrices_[chr_id].size()):
                matrices[chr_id].append(fm_to_array(matrices_[chr_id][k]))

    return matrices
