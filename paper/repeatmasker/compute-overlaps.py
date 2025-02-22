import tqdm
import numpy as np
import pandas as pd
import pyranges


bed = pyranges.read_bed('hg38_2020_rmsk.bed')


for cancer_type in ['colorectal', 'breast']:

    df = pd.read_csv(f'../../regions/metdecode/dmrs-methmarkerdb-{cancer_type}.csv', header=None, names=['chr', 'start', 'end', 'cat'])

    overlaps = []
    for chr_id, start, end in tqdm.tqdm(list(zip(df['chr'], df['start'], df['end']))):
        if chr_id < 23:
            chr_name = f'chr{chr_id}'
        elif chr_id == 23:
            chr_name = 'chrX'
        elif chr_id == 24:
            chr_name = 'chrY'
        else:
            chr_name = str(chr_id)
        start, end = int(start), int(end)
        genomic_interval = pyranges.from_dict({
            'Chromosome': [chr_name],
            'Start': [start],
            'End': [end]
        })

        ov = bed.overlap(genomic_interval)
        x = np.zeros(end - start + 1, dtype=bool)
        if len(ov.df) > 0:
            for _start, _end in zip(ov.df['Start'], ov.df['End']):
                x[max(0, _start - start):min(len(x), _end - end)] = True
        overlap = int(np.sum(x))
        print(overlap, end - start + 1, overlap / float(end - start + 1))
        assert overlap <= end - start + 1
        overlaps.append(overlap)

    df['repeat-masker-overlap'] = np.asarray(overlaps)

    df.to_csv(f'../../regions/metdecode/dmrs-methmarkerdb-{cancer_type}-with-rmo.csv', sep=',', header=False, index=False)
