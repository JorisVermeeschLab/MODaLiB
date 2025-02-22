from typing import Tuple
import os
import re
import argparse
import json
import tqdm

import numpy as np

from modalibpy.mbias import compute_bounds


folders = [
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20230823_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20230928_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20231004_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20231019_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20231107_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20231120_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20231130_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20231205_Angelica/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240131_Angelica/methy',
    '/lustre1/project/stg_00019/research/qfu/Project/methy_data_stg19_tmp/20240206/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/202301/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/202306/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/202307/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/202308/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20230928/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20230928_rerun_fix_contamination/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240111/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240119/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240212/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240212.part2/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240214.pt3/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240214.pt4/methy',
    '/lustre1/project/stg_00124/methylation_data/preprocessing/20240228/methy',
    '/lustre1/project/stg_00019/dhanya_research/emseq/apr2021/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/aug2021/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/dec2020/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/feb2022/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/feb2023/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/jan2022/align/mbias_image',
    '/lustre1/project/stg_00019/dhanya_research/emseq/jan2023/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/july2021/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/mar2023/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/march2022/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/may2021/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/nov2020/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/nov2021/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/nov2022/align',
    '/lustre1/project/stg_00019/dhanya_research/emseq/reanalysis/align'
]


REGEX = re.compile(r'(GC[0-9]+)(_[ACGT]+-[ACGT]+)?(_.+)?.mbias_(OB|OT).svg')


results = {}

for folder in folders:

    for filename in tqdm.tqdm(os.listdir(folder), desc=folder):

        filepath = os.path.join(folder, filename)
        if not (filepath.endswith('.mbias_OT.svg') or filepath.endswith('.mbias_OB.svg')):
            continue

        result = REGEX.match(filename)
        result = result.groups()
        gc_code = result[0]
        if result[1] is not None:
            gc_code = gc_code + result[1]
        strand = result[3]

        if gc_code not in results:
            results[gc_code] = {}
        curves, bounds = compute_bounds(filepath)
        if strand == 'OT':
            results[gc_code]['OT-bounds'] = bounds
            results[gc_code]['OT-R1'] = curves[0]
            results[gc_code]['OT-R2'] = curves[1]
        else:
            results[gc_code]['OB-bounds'] = bounds
            results[gc_code]['OB-R1'] = curves[0]
            results[gc_code]['OB-R2'] = curves[1]

    with open('inclusion-bounds.tsv', 'w') as f:
        for gc_code in results.keys():
            ot_str = ','.join([str(x) for x in results[gc_code]['OT-bounds']])
            ob_str = ','.join([str(x) for x in results[gc_code]['OB-bounds']])
            f.write(f'{gc_code}\t{ot_str}\t{ob_str}\n')

    with open('mbias-curves.json', 'w') as f:
        json.dump(results, f)
