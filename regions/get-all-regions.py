
regions = []

chr_id = 0
start = 1
region_id = 1
with open('map_hg38_50kb.wig', 'r') as f:
    for line in f.readlines():
        line = line.rstrip()
        elements = line.split()
        if len(elements) > 1:
            assert elements[1].startswith('chrom')
            chr_name = elements[1].split('=')[1]
            if chr_name.startswith('chr'):
                chr_name = chr_name[3:]
            if chr_name == 'X':
                chr_id = 23
            elif chr_name == 'Y':
                chr_id = 24
            else:
                chr_id = int(chr_name)
            start = 1
        else:
            end = start + 49999
            mappability = float(elements[0])
            assert 0 <= mappability <= 1
            if mappability >= 0.8:
                regions.append([chr_id, start, end, region_id])
            start += 50000
            region_id += 1

with open('all-regions.csv', 'w') as f:
    for i, (chr_id, start, end, region_id) in enumerate(regions):
        f.write(f'{chr_id},{start},{end},region{region_id}\n')
