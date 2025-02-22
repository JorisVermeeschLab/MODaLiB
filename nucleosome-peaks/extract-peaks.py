import pyBigWig
import tqdm


# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378
bb = pyBigWig.open('GSE71378_CA01.bb')
chr_sizes = bb.chroms()

with open('GSE71378_CA01.peaks.txt', 'w') as f:
    for chr_id in list(range(1, 23)) + ['X', 'Y']:
        chr_size = chr_sizes[f'chr{chr_id}']
        batch_size = 100000
        for batch_start in tqdm.tqdm(range(1, chr_size + batch_size + 1, batch_size), desc=f'chr{chr_id}'):
            batch_end = min(batch_start + batch_size, chr_size)
            if batch_end < batch_start:
                continue

            try:
                entries = bb.entries(f'chr{chr_id}', batch_start, batch_end)
            except RuntimeError as e:
                print(f'chr{chr_id}', batch_start, batch_end, chr_size)
                raise e

            if entries is not None:
                for entry in entries:
                    start, end, _ = entry
                    pos = (int(end) + int(start)) // 2
                    f.write(f'chr{chr_id}\t{pos}\n')
