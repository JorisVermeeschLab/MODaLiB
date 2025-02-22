from typing import List, Tuple


def convert_wig_file(filepath: str, out_filepath: str) -> None:

    regions = []

    chr_id = 0
    start = 1
    region_id = 1
    with open(filepath, 'r') as f:
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
                value = float(elements[0])
                regions.append([chr_id, start, end, value])
                start += 50000
                region_id += 1

    with open(out_filepath, 'w') as f:
        for chr_id, start, end, value in regions:
            f.write(f'{chr_id}\t{start}\t{end}\t{value}\n')


convert_wig_file('gc_hg38_50kb.wig', 'gc_hg38_50kb.bed')
