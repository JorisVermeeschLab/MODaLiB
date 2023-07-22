# -*- coding: utf-8 -*-
#
# core.pxd
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

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp cimport bool

import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef extern from "modalib/bounds.h" namespace "modalib":

    cdef int coords_to_partition(cnp.uint8_t chr_id, long pos, size_t n_partitions)


cdef extern from "modalib/region_list.h" namespace "modalib":

    cdef cppclass _RegionList:
        _RegionList()
        void load(const string &filepath)


cdef extern from "modalib/read_info.h" namespace "modalib":

    cdef struct ReadInfo:
        cnp.uint8_t start_chr_id
        cnp.uint32_t start_pos
        cnp.uint16_t n_methylated
        cnp.uint16_t n_unmethylated
        cnp.uint16_t gc_count
        cnp.uint16_t s_count
        cnp.uint16_t size
        float nuc_score


cdef extern from "modalib/sample.h" namespace "modalib":

    cdef cppclass _Sample:
        _Sample() nogil
        _Sample(string filepath) nogil
        void operator=(const _Sample&) nogil
        ReadInfo& get() nogil
        void next() nogil
        bool eof() nogil


cdef extern from "modalib/features.h" namespace "modalib":
    
    cdef cppclass FeatureMatrix:
        FeatureMatrix() nogil
        FeatureMatrix(size_t n, size_t m) nogil
        FeatureMatrix(FeatureMatrix& other) nogil
        double operator()(int i, int j) nogil
        size_t n_rows() nogil
        size_t n_cols() nogil

    cdef unordered_map[cnp.uint8_t, vector[FeatureMatrix]] extract_all_features(
            vector[_Sample]& samples,
            vector[bool] chr_mask,
            unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
    ) nogil


cdef extern from "modalib/kernel.h" namespace "modalib":

    cdef cppclass KernelMatrix:
        KernelMatrix() nogil
        KernelMatrix(size_t n, size_t m) nogil
        KernelMatrix(KernelMatrix& other) nogil
        double operator()(int i, int j) nogil
        size_t n_rows() nogil
        size_t n_cols() nogil

    cdef unordered_map[cnp.uint8_t, vector[KernelMatrix]] compute_kernel_matrices(
            vector[_Sample]& samples,
            vector[bool]& chr_mask,
            unordered_map[cnp.uint8_t, vector[long]]& nuc_pos
    ) nogil
