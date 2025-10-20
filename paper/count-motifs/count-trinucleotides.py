from typing import Optional, List, Dict

import os
import sys
import tqdm
import random
import time
import json

import numpy as np
import pysam

from modalib.core import Sample, ReferenceNucMap, chr_id_to_chr_name


ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, '..', 'step1-bam-processing', 'out_')
OUT_DIR = os.path.join(ROOT, 'outputs')
REF_GENOME_PATH = os.path.join(ROOT, '..', 'references', 'hg38.fa')
nuc_filepath = os.path.join(ROOT, '..', 'nucleosome-peaks', 'GSE71378_CA01.peaks.hg38.txt')
os.makedirs(OUT_DIR, exist_ok=True)


ref = pysam.FastaFile(REF_GENOME_PATH)
TRINUCLEOTIDES = {
    'AAA': 0, 'AAC': 1, 'AAG': 2, 'AAT': 3, 'ACA': 4, 'ACC': 5, 'ACG': 6, 'ACT': 7, 'AGA': 8, 'AGC': 9, 'AGG': 10, 'AGT': 11, 'ATA': 12, 'ATC': 13, 'ATG': 14, 'ATT': 15,
    'CAA': 16, 'CAC': 17, 'CAG': 18, 'CAT': 19, 'CCA': 20, 'CCC': 21, 'CCG': 22, 'CCT': 23, 'CGA': 24, 'CGC': 25, 'CGG': 26, 'CGT': 27, 'CTA': 28, 'CTC': 29, 'CTG': 30, 'CTT': 31,
    'GAA': 32, 'GAC': 33, 'GAG': 34, 'GAT': 35, 'GCA': 36, 'GCC': 37, 'GCG': 38, 'GCT': 39, 'GGA': 40, 'GGC': 41, 'GGG': 42, 'GGT': 43, 'GTA': 44, 'GTC': 45, 'GTG': 46, 'GTT': 47,
    'TAA': 48, 'TAC': 49, 'TAG': 50, 'TAT': 51, 'TCA': 52, 'TCC': 53, 'TCG': 54, 'TCT': 55, 'TGA': 56, 'TGC': 57, 'TGG': 58, 'TGT': 59, 'TTA': 60, 'TTC': 61, 'TTG': 62, 'TTT': 63
}

nucmap = ReferenceNucMap(nuc_filepath, distance=True)


def add_trinucleotide_frequencies(freqs: np.ndarray, sequence: str) -> None:
    for i in range(len(sequence) - 2):
        tri = sequence[i:i+3]
        if tri in TRINUCLEOTIDES:
            freqs[TRINUCLEOTIDES[tri]] += 1


def count_trinucleotides(filepath, timeout: float = 5 * 60) -> np.ndarray:
    sample = Sample(filepath)
    n = len(sample)
    freqs = np.zeros(64, dtype=np.uint64)
    t0 = time.time()
    n_frags = 0
    while True:
        i = random.randint(0, n - 1)
        frag = sample[i]

        chr_id = frag['start_chr_id']
        chr_name = chr_id_to_chr_name(chr_id)
        if chr_name in {'chr0', 'EBV'}:
            continue
        sequence = ref.fetch(chr_name, frag['start_pos'] - 1, frag['end_pos']).upper()

        add_trinucleotide_frequencies(freqs, sequence)
        n_frags += 1

        if time.time() - t0 >= timeout:
            break
    print(n_frags)

    return freqs


def count_nucleotides_per_position(filepath: str, timeout: float = 3600) -> Dict[str, np.ndarray]:
    sample = Sample(filepath)
    n = len(sample)
    rel_freqs_5p = np.zeros((64, 101), dtype=np.uint64)
    abs_freqs_5p = np.zeros((64, 500), dtype=np.uint64)
    rel_freqs_3p = np.zeros((64, 101), dtype=np.uint64)
    abs_freqs_3p = np.zeros((64, 500), dtype=np.uint64)
    nuc_abs_freqs = np.zeros((64, 500), dtype=np.uint64)
    t0 = time.time()
    n_frags = 0

    i = 0
    while True:
        #i = random.randint(0, n - 1)
        frag = sample[i]
        i += 1

        chr_id = frag['start_chr_id']
        chr_name = chr_id_to_chr_name(chr_id)
        if chr_name in {'chr0', 'EBV'}:
            continue
        sequence = ref.fetch(chr_name, frag['start_pos'] - 1, frag['end_pos']).upper()

        pos = frag['start_pos']

        for i in range(len(sequence)):
            tri = sequence[i:i+3]
            if tri not in TRINUCLEOTIDES:
                continue
            t = TRINUCLEOTIDES[tri]

            score = int(nucmap.score(chr_id, pos + i))
            if 0 <= score < 500:
                nuc_abs_freqs[t, score] += 1

                j = int(100 * i / len(sequence))
                if i < 500:
                    abs_freqs_5p[t, i] += 1
                rel_freqs_5p[t, j] += 1
                i_inv = len(sequence) - i - 1
                j_inv = int(100 * i_inv / len(sequence))
                if i_inv < 500:
                    abs_freqs_3p[t, i_inv] += 1
                rel_freqs_3p[t, j_inv] += 1

        print(n_frags, flush=True)
        n_frags += 1

        if time.time() - t0 >= timeout:
            break

    return {
        'rel_freqs_5p': rel_freqs_5p,
        'rel_freqs_3p': rel_freqs_3p,
        'abs_freqs_5p': abs_freqs_5p,
        'abs_freqs_3p': abs_freqs_3p,
        'nuc_abs_freqs': nuc_abs_freqs
    }


gc_code = 'GC131355'
filepath = os.path.join(DATA_DIR, f'{gc_code}.modalib')
np.savez(f'{gc_code}-perposition.npz', **count_nucleotides_per_position(filepath))
