CHR_LENGTHS = {
    "chr1" : (0, 4945),
    "chr2" : (1, 4860),
    "chr3" : (2, 3991),
    "chr4" : (3, 3826),
    "chr5" : (4, 3618),
    "chr6" : (5, 3418),
    "chr7" : (6, 3177),
    "chr8" : (7, 2926),
    "chr9" : (8, 2806),
    "chr10": (9, 2708),
    "chr11": (10, 2690),
    "chr12": (11, 2647),
    "chr13": (12, 2283),
    "chr14": (13, 2128),
    "chr15": (14, 2007),
    "chr16": (15, 1777),
    "chr17": (16, 1576),
    "chr18": (17, 1523),
    "chr19": (18, 1277),
    "chr20": (19, 1249),
    "chr21": (20, 939),
    "chr22": (21, 994),
    "chrX" : (22, 3099),
    "chrY" : (23, 1156),
    "chrM" : (24, 1),
    "EBV"  : (25, 4)
}

with open('bins-gip.txt', 'w') as f:
    for key in CHR_LENGTHS.keys():
        n = CHR_LENGTHS[key][1]
        for i in range(n):
            f.write(f'{key}, {i*50000+1}, {(i+1)*50000+1}\n')


with open('ncbi_dataset.tsv', 'r') as f:
    lines = f.readlines()[1:]
regions = []
processed = set()
for line in lines:
    line = line.rstrip()
    elements = line.split('\t')
    if len(elements) == 10:
        if elements[0] in processed:
            continue
        processed.add(elements[0])
        ncbi = elements[0]
        gene_id = elements[1]
        chr_name = elements[6]
        start = int(elements[8])
        end = int(elements[9])

        if chr_name == 'X,Y':
            continue
        elif chr_name == 'X':
            idx = 23
            chr_name = 'chrX'
        elif chr_name == 'Y':
            idx = 24
            chr_name = 'chrY'
        elif chr_name == 'MT':
            idx = 25
            chr_name = 'chrM'
        elif chr_name == 'Un':
            continue
        else:
            idx = f'chr{int(chr_name)}'
        idx = idx * 1e+10 + start

        regions.append((idx, chr_name, start, end, ncbi, gene_id))
regions = list(sorted(regions))
with open('genes.txt', 'w') as f:
    for _, chr_name, start, end, ncbi, gene_id in regions:
        f.write(f'{chr_name}, {start}, {end}, {ncbi}, {gene_id}\n')
