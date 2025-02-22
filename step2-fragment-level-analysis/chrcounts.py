import re
import os
import time
import argparse
import tqdm
import collections

import numpy as np
from modalib.core import Sample


ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(ROOT, '../step1-bam-processing/out')

GC_CODE_PATTERN = re.compile(r'^GC([0-9]{6,})(-([A-Z0-9]{5}))?')


parser = argparse.ArgumentParser()
parser.add_argument(
    '-gc_code',
    type=str,
    default='all',
    help='GC code'
)

args = parser.parse_args()
gc_code = args.gc_code


def process_file(filepath: str) -> None:
    sample = Sample(filepath)
    counter = collections.Counter()
    for i, fragment in enumerate(tqdm.tqdm(sample)):
        counter[fragment['start_chr_id']] += 1
    print(counter)


def gc_code_from_filename(filename: str) -> str:
    return re.match(GC_CODE_PATTERN, filename).group(0)


if gc_code == 'all':
    for filename in os.listdir(DATA_PATH):
        if not filename.endswith('.modalib'):
            continue
        try:
            gc_code = gc_code_from_filename(filename)
        except AttributeError:
            continue
        sample = Sample(os.path.join(DATA_PATH, filename))
        n = len(sample)
        print(sample[n - 2]['start_chr_id'], gc_code)
else:
    sample = Sample(os.path.join(DATA_PATH, f'{gc_code}.modalib'))
    for frag in sample:
        print(frag)
        break
        print(sample[0])
    n = len(sample)
    print(sample[n - 2])
