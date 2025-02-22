from typing import Optional, Tuple, List, Dict
import re
import os
import time
import json
import random
import logging

import numpy as np

import pysam

from modalib.core import Sample, pairwise_kernels, extract_features, compute_mss_matrix, compute_mss_markers_matrix


def process_bam_file(gc_code: str, filepath: str, override_: bool = False) -> bool:
    out_filepath = os.path.join(OUT_DIR, f'{gc_code}.modalib')
    if override_ or (not os.path.exists(out_filepath)):
        if gc_code not in bounds:
            gc_code = gc_code.split('_')[0]
        if gc_code not in bounds:
            print(f'Unknown inclusion bounds for {gc_code}')
            ref_genome, _ = get_reference(filepath)
            cmd = f'{EXE_PATH} {ref_genome} {filepath} -o {out_filepath} --cconv'
            try:
                os.system(cmd)
            except Exception as e:
                print(e)
                logging.exception(e)
                return False
        else:
            ot, ob = bounds[gc_code]
            print(f'Processing {gc_code}')
            ref_genome, _ = get_reference(filepath)
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
        gc_code: str,
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

        if gc_code not in mbias_curves:
            gc_code = gc_code.split('_')[0]
        if gc_code in mbias_curves:
            for key in ['OT-R1', 'OT-R2', 'OB-R1', 'OB-R2']:
                matrices[key] = np.asarray(mbias_curves[gc_code][key], dtype=float)

        np.savez(out_filepath, **matrices)


for gc_code in GC_CODES:
    if (gc_code not in filepath_dict) and ('_' in gc_code):
        gc_code = gc_code.split('_')[0]

    if gc_code in os.listdir(MUTEX_DIR):
        continue
    with open(os.path.join(MUTEX_DIR, gc_code), 'w') as f:
        f.write('-')

    if gc_code not in filepath_dict:
        continue  # TODO

    filepath = filepath_dict[gc_code]
    out_folder = os.path.join(GC_PARAGON_LOCATION, 'out')
    os.makedirs(out_folder, exist_ok=True)

    gc_gc_code = os.path.basename(filepath).replace('.bam', '')
    gc_filepath = os.path.join(out_folder, gc_gc_code, f'{gc_gc_code}.GCtagged.bam')

    if not os.path.exists(gc_filepath):
        _, ref_genome = get_reference(filepath)

        cmd = f'python {GC_PARAGON_SCRIPT} --bam {filepath} -o {out_folder} --output-bam --two-bit-reference-genome {ref_genome}'
        try:
            os.system(cmd)
        except Exception as e:
            print(e)
            logging.exception(e)
            continue

    modalib_filepath = os.path.join(OUT_DIR, f'{gc_code}.modalib')
    if not os.path.exists(modalib_filepath):
        succeeded = process_bam_file(gc_code, gc_filepath)
        if not succeeded:
            continue


    #for prefix, kwargs in zip(['', 'downsampled-', 'nogc-', 'nogc-downsampled-'], [{}, {'target_n_fragments': 4000000}, {'weighting': False}, {'weighting': False, 'target_n_fragments': 4000000}]):
    #for prefix, kwargs in zip(['downsampled-'], [{'target_n_fragments': 4000000}]):
    #for prefix, kwargs in zip(['', 'downsampled-'], [{}, {'target_n_fragments': 4000000}]):
    for prefix, kwargs in zip([''], [{}]):
        
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}methmarkerdbbreast-noXY', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'metdecode', 'dmrs-methmarkerdb-breast.csv'),
            **kwargs
        )

        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}ctcf-noXY', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'ctcf-no-XY.csv'),
            **kwargs
        )
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}50kb-noXY', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'all-regions-no-XY.csv'),
            **kwargs
        )
        for key in ['promoters', 'tss', 'genebodies']:
            process_modalib_file(
                modalib_filepath,
                os.path.join(OUT_DIR_2, f'{prefix}{key}-noXY', f'{gc_code}.npz'),
                gc_code,
                os.path.join(ROOT, 'regions', f'{key}-no-XY.csv'),
                **kwargs
            )
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}metdecode-noXY', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'metdecode', 'dmrs.csv'),
            **kwargs
        )
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}methmarkerdbcolorectal-noXY', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'metdecode', 'dmrs-methmarkerdb-colorectal.csv'),
            **kwargs
        )
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}umac', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'umac.csv'),
            **kwargs
        )
        for modality in ['cin', 'fsi', 'mfr']:
            process_modalib_file(
                modalib_filepath,
                os.path.join(OUT_DIR_2, f'{prefix}themis-{modality}', f'{gc_code}.npz'),
                gc_code,
                os.path.join(ROOT, 'regions', 'themis', f'themis-{modality}-regions.csv'),
                **kwargs
            )

        """
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}mtintegrity', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'mt-integrity.csv'),
            **kwargs
        )
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}mt', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'mt.csv'),
            **kwargs
        )
        process_modalib_file(
            modalib_filepath,
            os.path.join(OUT_DIR_2, f'{prefix}nucleases', f'{gc_code}.npz'),
            gc_code,
            os.path.join(ROOT, 'regions', 'nucleases.csv'),
            **kwargs
        )
        """

    if os.path.exists(os.path.join(MUTEX_DIR, gc_code)):
        os.remove(os.path.join(MUTEX_DIR, gc_code))


print(f'Finished.')
