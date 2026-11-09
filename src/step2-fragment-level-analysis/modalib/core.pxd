# -*- coding: utf-8 -*-
#
# core.pxd
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=True

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.map cimport map
from libcpp cimport bool
from libc.time cimport time_t

import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef extern from "modalib/utils.h" namespace "modalib":
    
    ctypedef long long genome_size_t
    genome_size_t combine_coords(cnp.uint8_t chr_id, cnp.uint32_t pos)


cdef extern from "modalib/bounds.h" namespace "modalib":

    cdef int coords_to_partition(cnp.uint8_t chr_id, long pos, size_t n_partitions)


cdef extern from "modalib/region_list.h" namespace "modalib":

    cdef cppclass _RegionList:
        _RegionList()
        void load(const string &filepath)


cdef extern from "modalib/fragment_info.h" namespace "modalib":

    cdef struct FragmentInfo:
        cnp.uint8_t chr_id
        cnp.uint32_t start_pos
        cnp.uint32_t end_pos
        cnp.uint16_t size
        cnp.uint16_t insert_size
        cnp.uint8_t n_methylated_cpg_r1
        cnp.uint8_t n_unmethylated_cpg_r1
        cnp.uint8_t n_methylated_cpg_r2
        cnp.uint8_t n_unmethylated_cpg_r2
        cnp.uint8_t n_methylated_cpg
        cnp.uint8_t n_unmethylated_cpg
        cnp.uint8_t n_methylated_ch
        cnp.uint8_t n_unmethylated_ch
        cnp.uint16_t gc_count
        cnp.uint16_t s_count
        cnp.uint8_t end_motif_5p
        cnp.uint8_t end_motif_3p
        cnp.uint8_t strand
        cnp.uint8_t mapping_quality
        float weight


cdef extern from "modalib/sample.h" namespace "modalib":

    cdef cppclass _Sample:
        _Sample() nogil
        _Sample(string filepath) nogil
        void operator=(const _Sample&) nogil
        FragmentInfo& get() nogil
        void go_to_chunk(long i) nogil
        void next() nogil
        bool eof() nogil
        long get_length() nogil


cdef extern from "modalib/nucleosomes.h" namespace "modalib":

    cdef cppclass RefNucMap:
        RefNucMap() nogil
        RefNucMap(const string filepath, bool distance) nogil
        RefNucMap(const string filepath) nogil
        float phase_shift_score_wrt_reference(cnp.uint8_t chr_id, unsigned long x) nogil


cdef extern from "modalib/features.h" namespace "modalib":
    
    cdef cppclass FeatureMatrix:
        FeatureMatrix() nogil
        FeatureMatrix(size_t n, size_t m) nogil
        FeatureMatrix(size_t n, size_t m, bool) nogil
        FeatureMatrix(FeatureMatrix& other) nogil
        double operator()(int i, int j) nogil
        size_t n_rows() nogil
        size_t n_cols() nogil

    cdef vector[FeatureMatrix] extract_all_features(
            vector[_Sample]& samples,
            RefNucMap& nuc_map,
            _RegionList& regions,
            double target_n_fragments,
            bool weighting
    ) nogil


cdef extern from "modalib/microsatellites.h" namespace "modalib":

    cdef cppclass MSS:
        string filepath
        map[genome_size_t, unordered_map[int, size_t]] data
        map[genome_size_t, size_t] counts
        map[genome_size_t, string] motifs
        MSS()
        MSS(string filepath)
        void extract_features(int i, FeatureMatrix& features, FeatureMatrix& counts)
        void extract_markers(int i, FeatureMatrix& features, FeatureMatrix& counts, const unordered_map[genome_size_t, size_t]& markers_pos)
        void close()


cdef extern from "modalib/dms.h" namespace "modalib":

    cdef void find_differentially_methylated_sites(
            const vector[string]& filepaths,
            const vector[float]& y,
            const string& out_filepath)


cdef extern from "modalib/kernel.h" namespace "modalib":

    cdef cppclass KernelMatrix:
        KernelMatrix() nogil
        KernelMatrix(size_t n) nogil
        KernelMatrix(KernelMatrix& other) nogil
        double operator()(int i, int j) nogil
        size_t n_rows() nogil
        size_t n_cols() nogil

    cdef vector[KernelMatrix] compute_kernel_matrices(
            vector[_Sample]& samples,
            RefNucMap& nuc_map,
            _RegionList& regions
    ) nogil
