from typing import Optional, Tuple, List, Dict
import re
import os
import time
import json
import random
import logging

import numpy as np

import pysam

from modalib.core import Sample, extract_features


ROOT = os.path.dirname(os.path.abspath(__file__))
EXE_PATH = os.path.join(ROOT, 'step1-bam-processing', 'BAMConvertor')

OUT_DIR = os.path.join(ROOT, 'step1-bam-processing', 'out_')
os.makedirs(OUT_DIR, exist_ok=True)
OUT_DIR_2 = os.path.join(ROOT, 'outputs')
os.makedirs(OUT_DIR_2, exist_ok=True)


nuc_filepath = "/lustre1/project/stg_00019/research/Antoine/MODaLiB/nucleosome-peaks/GSE71378_CA01.peaks.hg38.txt"

filepath = "/lustre1/project/stg_00149/MICADO/methylation_data/preprocessing/20240214.pt3/bam/GC146766_AACCAGAG-CTTCGGTT.dupmark.merged.bam"


ref_genome = "/staging/leuven/stg_00149/MICADO/methylation_data/resource/reference/hg38_lambda_puc19/GCA_000001405.15_GRCh38_no_alt_analysis_set_lambda_puc19_combined.fa"
ref_genome_2bit = "/lustre1/project/stg_00019/research/Antoine/MODaLiB/references/hg38-lambda-puc19.2bit"

GC_PARAGON_SCRIPT = '/lustre1/project/stg_00019/research/Antoine/GCparagon/src/GCparagon/correct_GC_bias.py'


def process_bam_file(filepath: str, out_filepath: str, override_: bool = False, bounds: Tuple[str, str] = ("0,0,0,0", "0,0,0,0")) -> bool:
    if override_ or (not os.path.exists(out_filepath)):
        ot, ob = bounds
        cmd = f'{EXE_PATH} {ref_genome} {filepath} -o {out_filepath} --cconv --OT {ot} --OB {ob}'
        try:
            os.system(cmd)
        except Exception as e:
            print(e)
            logging.exception(e)
            return False
    return True


def process_modalib_file(
        in_filepath: str,
        out_filepath: str,
        region_list_filepath: str,
        target_n_fragments: Optional[int] = None,
        override_: bool = False,
        lightweight: bool = False,
        weighting: bool = True
    ) -> None:

    if not os.path.exists(in_filepath):
        return

    os.makedirs(os.path.dirname(out_filepath), exist_ok=True)

    if (not os.path.exists(out_filepath)) or override_:
        if target_n_fragments is not None:
            matrices = extract_features([in_filepath], nuc_filepath, region_list_filepath, target_n_fragments=target_n_fragments, weighting=weighting)
        else:
            matrices = extract_features([in_filepath], nuc_filepath, region_list_filepath, weighting=weighting)

        for key in matrices:
            matrices[key] = np.squeeze(matrices[key])

        if lightweight:
            for key in list(matrices.keys()):
                if key not in {'Fragment length counts', 'End-motif probability (5 prime)', 'End-motif probability (3 prime)', 'Nucleosome score histogram', 'R1 methylation', 'R2 methylation'}:
                    del matrices[key]

        np.savez(out_filepath, **matrices)


filename = os.path.splitext(os.path.basename(filepath))[0]

out_folder = os.path.dirname(os.path.abspath(GC_PARAGON_SCRIPT))
os.makedirs(out_folder, exist_ok=True)
gc_filepath = os.path.join(out_folder, f'{filename}.GCtagged.bam')

if not os.path.exists(gc_filepath):

    cmd = f'python {GC_PARAGON_SCRIPT} --bam {filepath} -o {out_folder} --output-bam --two-bit-reference-genome {ref_genome_2bit}'
    os.system(cmd)

modalib_filepath = os.path.join(OUT_DIR, f'{filename}.modalib')
if not os.path.exists(modalib_filepath):
    process_bam_file(gc_filepath, modalib_filepath)

for prefix, kwargs in zip([''], [{}]):

    process_modalib_file(
        modalib_filepath,
        os.path.join(OUT_DIR_2, f'{prefix}methmarkerdbbreast-noXY', f'{filename}.npz'),
        os.path.join(ROOT, 'regions', 'metdecode', 'dmrs-methmarkerdb-breast.csv'),
        **kwargs
    )

    process_modalib_file(
        modalib_filepath,
        os.path.join(OUT_DIR_2, f'{prefix}ctcf-noXY', f'{filename}.npz'),
        os.path.join(ROOT, 'regions', 'ctcf-no-XY.csv'),
        **kwargs
    )
    process_modalib_file(
        modalib_filepath,
        os.path.join(OUT_DIR_2, f'{prefix}50kb-noXY', f'{filename}.npz'),
        os.path.join(ROOT, 'regions', 'all-regions-no-XY.csv'),
        **kwargs
    )
    for key in ['promoters', 'tss', 'genebodies']:
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}{key}-noXY', f'{filename}.npz'),
            os.path.join(ROOT, 'regions', f'{key}-no-XY.csv'),
            **kwargs
        )
    process_modalib_file(
        modalib_filepath,
        os.path.join(OUT_DIR_2, f'{prefix}metdecode-noXY', f'{filename}.npz'),
        os.path.join(ROOT, 'regions', 'metdecode', 'dmrs.csv'),
        **kwargs
    )
    process_modalib_file(
        modalib_filepath,
        os.path.join(OUT_DIR_2, f'{prefix}methmarkerdbcolorectal-noXY', f'{filename}.npz'),
        os.path.join(ROOT, 'regions', 'metdecode', 'dmrs-methmarkerdb-colorectal.csv'),
        **kwargs
    )


print(f'Finished.')
