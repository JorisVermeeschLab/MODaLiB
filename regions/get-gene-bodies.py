import tqdm
import gzip
import numpy as np
import pandas as pd


# https://github.com/shendurelab/cfDNA/tree/master/data/Ensembl_v75
transcript_ids = set()
with gzip.open('ensemblv75_canonicalTranscriptIDs.protein_coding.lst.gz', 'r') as f:
    for line in tqdm.tqdm(f, desc='Reading .lst.gz file'):
        line = line.decode('utf-8').strip()
        transcript_ids.add(line)

# https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/
gene_bodies = {}
with gzip.open('Homo_sapiens.GRCh38.110.gtf.gz', 'r') as f:

    for line in tqdm.tqdm(f, desc='Reading .gtf.gz file'):
        line = line.decode('utf-8').rstrip()

        if line.startswith('#'):
            continue

        elements = line.split('\t')
        fields = elements[8].split(';')

        transcript_id = ''
        gene_id = ''
        gene_name = ''
        for field in fields:
            if 'transcript_id' in field:
                transcript_id = field.replace('transcript_id', '').replace('"', '').strip()
            elif 'gene_id' in field:
                gene_id = field.replace('gene_id', '').replace('"', '').strip()
            elif 'gene_name' in field:
                gene_name = field.replace('gene_name', '').replace('"', '').strip()

        if (elements[2] == 'transcript') and (transcript_id in transcript_ids):

            if elements[0] == 'X':
                chr_id = 23
            elif elements[0] == 'Y':
                chr_id = 24
            elif elements[0] == 'MT':
                chr_id = 25
            elif elements[0].isnumeric():
                chr_id = int(elements[0])
            elif elements[0].startswith('K'):
                continue
            else:
                raise ValueError(elements[0])

            start = int(elements[3])
            end = int(elements[4])
            coord = 300000000 * chr_id + start

            # Keep only first 10kb of gene body
            end = min(end, start + 9999)

            gene_bodies[coord] = (chr_id, start, end, gene_name)

coords = np.asarray(list(gene_bodies.keys()), dtype=np.uint64)
coords = np.sort(coords)

with open('gene-bodies.csv', 'w') as f:
    for coord in coords:
        chr_id, start, end, gene_name = gene_bodies[int(coord)]
        f.write(f'{chr_id},{start},{end},{gene_name}\n')

with open('gene-bodies-no-XY.csv', 'w') as f:
    for coord in coords:
        chr_id, start, end, gene_name = gene_bodies[int(coord)]
        if chr_id <= 22:
            f.write(f'{chr_id},{start},{end},{gene_name}\n')

with open('tss.csv', 'w') as f:
    for coord in coords:
        chr_id, start, _, gene_name = gene_bodies[int(coord)]
        start = max(1, start - 50)
        end = start + 50
        f.write(f'{chr_id},{start},{end},{gene_name}\n')

with open('tss-no-XY.csv', 'w') as f:
    for coord in coords:
        chr_id, start, _, gene_name = gene_bodies[int(coord)]
        if chr_id <= 22:
            start = max(1, start - 50)
            end = start + 50
            f.write(f'{chr_id},{start},{end},{gene_name}\n')
