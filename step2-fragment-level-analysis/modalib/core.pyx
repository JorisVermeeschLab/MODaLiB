# -*- coding: utf-8 -*-
#
# core.pyx
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=True

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

def uint8_motif_to_str(motif: int) -> str:
    BASES = ['A', 'T', 'C', 'G']
    s = BASES[(motif & 0xC0) >> 6]
    s += BASES[(motif & 0x30) >> 4]
    s += BASES[(motif & 0x0C) >> 2]
    s += BASES[(motif & 0x03)]
    return s


def chr_id_to_chr_name(chr_id: int) -> str:
    if chr_id == 23:
        chr_name = 'chrX'
    elif chr_id == 24:
        chr_name = 'chrY'
    elif chr_id == 25:
        chr_name = 'chrM'
    elif chr_id == 26:
        chr_name = 'EBV'
    else:
        chr_name = f'chr{chr_id}'
    return chr_name


cdef class ReferenceNucMap:

    cdef string filepath
    cdef bool distance
    cdef RefNucMap base

    def __init__(self, filepath: str, distance: bool = False):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f'Could not open file {filepath}')
        self.filepath = <string>filepath.encode('utf-8')
        self.distance = distance
        self.base = RefNucMap(self.filepath, self.distance)

    def score(self, chr_id: int, pos: int) -> float:
        return self.base.phase_shift_score_wrt_reference(chr_id, pos)


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

    @staticmethod
    def frag_to_dict(FragmentInfo& frag_) -> dict:
        return dict(
            start_chr_id=frag_.chr_id,
            start_pos=frag_.start_pos,
            end_pos=frag_.end_pos,
            size=frag_.size,
            insert_size=frag_.insert_size,
            n_methylated_cpg_r1=frag_.n_methylated_cpg_r1,
            n_unmethylated_cpg_r1=frag_.n_unmethylated_cpg_r1,
            n_methylated_cpg_r2=frag_.n_methylated_cpg_r2,
            n_unmethylated_cpg_r2=frag_.n_unmethylated_cpg_r2,
            n_methylated_cpg=frag_.n_methylated_cpg,
            n_unmethylated_cpg=frag_.n_unmethylated_cpg,
            n_methylated_ch=frag_.n_methylated_ch,
            n_unmethylated_ch=frag_.n_unmethylated_ch,
            gc_count=frag_.gc_count,
            s_count=frag_.s_count,
            end_motif_5p_int=frag_.end_motif_5p,
            end_motif_3p_int=frag_.end_motif_3p,
            end_motif_5p=uint8_motif_to_str(frag_.end_motif_5p),
            end_motif_3p=uint8_motif_to_str(frag_.end_motif_3p),
            strand=frag_.strand,
            mapq=frag_.mapping_quality,
            weight=frag_.weight
        )

    def __iter__(self):
        cdef FragmentInfo frag_
        while not self.sample_.eof():
            frag_ = self.sample_.get()
            yield Sample.frag_to_dict(frag_)
            self.sample_.next()

    def __getitem__(self, i: int):
        self.sample_.go_to_chunk(i)
        return Sample.frag_to_dict(self.sample_.get())

    def get_name(self) -> str:
        return bytes(self.filepath.c_str()).decode('utf-8')

    def __len__(self) -> int:
        return self.sample_.get_length()


cdef cnp.ndarray to_array(KernelMatrix &km):
    cdef size_t n = km.n_rows()
    cdef size_t m = km.n_cols()
    cdef size_t i, j
    cdef cnp.float_t[:, :] arr = np.empty((n, m), dtype=float)
    with nogil:
        for i in range(n):
            for j in range(m):
                arr[i, j] = km(i, j)
    return np.asarray(arr)


cdef vector[_Sample] samples_from_filepaths(filepaths: List[str]):
    cdef vector[_Sample] samples
    for filepath in filepaths:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f'Could not locate {filepath}')
        samples.push_back(_Sample(filepath.encode('utf-8')))
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
        nuc_filepath: str,
        regions_filepath: str
) -> Dict[int, List[np.ndarray]]:

    # Check number of samples
    n = len(filepaths)
    if n <= 1:
        raise ValueError(f'Number of samples should be at least 2. Got {n}.')

    # Load reference nucleosome map
    if not os.path.exists(nuc_filepath):
        raise FileNotFoundError(f'Could not locate {nuc_filepath}')
    cdef RefNucMap nuc_map = RefNucMap(nuc_filepath.encode('utf-8'))

    # Load regions to be considered
    cdef _RegionList regions
    regions.load(regions_filepath.encode('utf-8'))

    # Compute kernel matrices
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[KernelMatrix] matrices_ = compute_kernel_matrices(samples, nuc_map, regions)

    # Convert output kernel matrices to dict of lists of numpy arrays
    matrices = []
    for k in range(matrices_.size()):
        matrices.append(to_array(matrices_.at(k)))

    return np.asarray(matrices)


# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


cdef cnp.ndarray fm_to_array(FeatureMatrix &km):
    cdef size_t n = km.n_rows()
    cdef size_t m = km.n_cols()
    cdef size_t i, j
    cdef cnp.float_t[:, :] arr = np.empty((n, m), dtype=float)
    with nogil:
        for i in range(n):
            for j in range(m):
                arr[i, j] = km(i, j)
    return np.asarray(arr)


def extract_features(
        filepaths: List[str],
        nuc_filepath: str,
        regions_filepath: str,
        target_n_fragments: float = 1e+25,
        weighting: bool = True
) -> Dict[str, np.ndarray]:

    # Check number of samples
    n = len(filepaths)
    if n < 1:
        raise ValueError(f'Number of samples should be at least 1. Got {n}.')

    # Load reference nucleosome map
    if not os.path.exists(nuc_filepath):
        raise FileNotFoundError(f'Could not find "{nuc_filepath}".')
    cdef RefNucMap nuc_map = RefNucMap(nuc_filepath.encode('utf-8'))

    # Load regions to be considered
    cdef _RegionList regions
    regions.load(regions_filepath.encode('utf-8'))

    # Compute feature matrices
    cdef vector[_Sample] samples = samples_from_filepaths(filepaths)
    cdef vector[FeatureMatrix] matrices_ = extract_all_features(samples, nuc_map, regions, target_n_fragments, weighting)

    # Convert output feature matrices to dict of list of numpy arrays
    matrices_dict = {
            'Nucleosome positioning': fm_to_array(matrices_.at(0)),
            'Fragment size index': fm_to_array(matrices_.at(1)),
            'Fragment length counts': fm_to_array(matrices_.at(2)),
            'CpG Methylation (first moment)': fm_to_array(matrices_.at(3)),
            'CpG Methylation (second moment)': fm_to_array(matrices_.at(4)),
            'CH Methylation (first moment)': fm_to_array(matrices_.at(5)),
            'Coverage': fm_to_array(matrices_.at(6)),
            'GC-content': fm_to_array(matrices_.at(7)),
            'CpG count': fm_to_array(matrices_.at(8)),
            'End-motif probability (5 prime)': fm_to_array(matrices_.at(9)),
            'End-motif probability (3 prime)': fm_to_array(matrices_.at(10)),
            'Mappability': fm_to_array(matrices_.at(11)),
            'Nucleosome score counts': fm_to_array(matrices_.at(12)),
            'Nucleosome score histogram': fm_to_array(matrices_.at(13)),
            'Fully methylated fragment ratio': fm_to_array(matrices_.at(14)),
            'CAFF': fm_to_array(matrices_.at(15)),
            'R1 methylation': fm_to_array(matrices_.at(16)),
            'R2 methylation': fm_to_array(matrices_.at(17)),
            'End-motif probability (THEMIS)': fm_to_array(matrices_.at(18)),
            'Proportion of long fragments': fm_to_array(matrices_.at(19)),
    }
    return matrices_dict


def compute_mss_matrix(filepaths: List[str]):
    cdef vector[string] filepaths_
    for filepath in filepaths:
        filepaths_.push_back(filepath.encode('utf-8'))

    cdef size_t n = filepaths_.size()
    cdef size_t m = 300
    cdef FeatureMatrix features = FeatureMatrix(n, m)
    cdef FeatureMatrix counts = FeatureMatrix(n, m, False)
    cdef MSS mss
    for i in range(n):
        mss = MSS(filepaths_.at(i))
        mss.extract_features(i, features, counts)
        mss.close()
    return fm_to_array(features), fm_to_array(counts)


def parse_chr_name(chr_name):
    if chr_name.startswith('chr'):
        chr_name = chr_name[3:]
    chr_name = chr_name.upper()
    if chr_name.isdigit():
        return int(chr_name)
    elif chr_name == 'X':
        return 23
    elif chr_name == 'Y':
        return 24
    elif chr_name == 'M':
        return 25
    elif chr_name == 'EBV':
        return 26
    else:
        raise NotImplementedError(f'Unknown chromosome name "{chr_name}"')


def compute_mss_markers_matrix(filepaths, marker_list_filepath):
    cdef vector[string] filepaths_
    for filepath in filepaths:
        filepaths_.push_back(filepath.encode('utf-8'))

    cdef cnp.uint8_t chr_id
    cdef genome_size_t combined_coords
    cdef unordered_map[genome_size_t, size_t] markers_pos
    with open(marker_list_filepath, 'r') as f:
        for line in f.readlines():
            elements = line.rstrip().split('\t')
            if len(elements) >= 2:
                chr_id = parse_chr_name(elements[0])
                pos = int(elements[1])
                combined_coords = combine_coords(chr_id, pos)
                markers_pos[combined_coords] = markers_pos.size()

    cdef size_t n = filepaths_.size()
    cdef size_t m = markers_pos.size()
    cdef FeatureMatrix features = FeatureMatrix(n, m)
    cdef FeatureMatrix counts = FeatureMatrix(n, m, False)
    cdef MSS mss
    for i in range(n):
        mss = MSS(filepaths_.at(i))
        mss.extract_markers(i, features, counts, markers_pos)
        mss.close()

    return fm_to_array(features), fm_to_array(counts)


def find_dms(filepaths: List[str], y: np.ndarray, out_filepath: str) -> None:
    cdef vector[string] filepaths_
    for filepath in filepaths:
        filepaths_.push_back(filepath.encode('utf-8'))

    cdef string out_filepath_ = out_filepath.encode('utf-8')

    cdef cnp.float_t[:] y_arr = y.astype(float)
    cdef vector[float] y_
    cdef int i
    with nogil:
        for i in range(y_arr.shape[0]):
            y_.push_back(y_arr[i])

    find_differentially_methylated_sites(filepaths_, y_, out_filepath_)
