import pandas as pd


df = pd.read_csv('MFR_matrix.txt', sep='\t')
regions = df['Window']
with open('themis-mfr-regions.csv', 'w') as f:
	for region in regions:
		left, right = region.split(':')
		chr_id = int(left.replace('chr', ''))
		start, end = right.split('-')
		start, end = int(start), int(end)
		f.write(f'{chr_id},{start},{end}\n')


import sys; sys.exit(0)


df = pd.read_csv('FSI_regions.txt', sep='\t')
regions = []
for chr_name, start, end in zip(df['seqnames'], df['start'], df['end']):
	chr_id = int(chr_name.replace('chr', ''))
	regions.append((chr_id, int(start), int(end)))
regions.sort(key=lambda x: 1e+9 * x[0] + x[1])
with open('themis-fsi-regions.csv', 'w') as f:
	for chr_id, start, end in regions:
		f.write(f'{chr_id},{start},{end}\n')


df = pd.read_csv('CIN_matrix.txt', sep='\t')
regions = df['Window']
with open('themis-cin-regions.csv', 'w') as f:
	for region in regions:
		left, right = region.split(':')
		chr_id = int(left.replace('chr', ''))
		start, end = right.split('-')
		start, end = int(start), int(end)
		f.write(f'{chr_id},{start},{end}\n')
