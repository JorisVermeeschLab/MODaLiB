import tqdm
import pandas as pd


df = pd.read_csv('CTCF_hg38_all.bed', sep='\t')

sites = []
for chr_name, start, end, p_value in tqdm.tqdm(zip(df['seqnames'], df['start'], df['end'], df['pvalue'])):
    if p_value < 1e-7:
        if chr_name.startswith('chr'):
            chr_name = chr_name[3:]
        if chr_name.isdigit():
            chr_id = int(chr_name)
        elif chr_name == 'X':
            chr_id = 23
        elif chr_name == 'Y':
            chr_id = 24
        else:
            raise NotImplementedError(f'Unknown contig "{chr_id}"')
        sites.append((chr_id, start, end))

with open('ctcf.csv', 'w') as f:
    for k, (chr_id, start, end) in enumerate(sites):
        f.write(f'{chr_id},{start},{end},{k}\n')

with open('ctcf-no-XY.csv', 'w') as f:
    for k, (chr_id, start, end) in enumerate(sites):
        if chr_id <= 22:
            f.write(f'{chr_id},{start},{end},{k}\n')
