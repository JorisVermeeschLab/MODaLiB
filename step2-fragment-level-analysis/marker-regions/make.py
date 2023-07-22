import os


class Region:

    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.chr_id = self.parse_chr_id(self.chr_name)
        self.start = int(start)
        self.end = int(end)

    def parse_chr_id(self, chr_name):
        if chr_name.startswith('chr'):
            chr_name = chr_name[3:]
        if len(chr_name) == 0:
            return 0
        if chr_name.isdigit():
            return int(chr_name)
        if chr_name[0] == 'X':
            return 23
        if chr_name[0] == 'Y':
            return 24
        if chr_name[0] == 'M':
            return 25
        return 0

    def __cmp__(self, other):
        if self.chr_id < other.chr_id:
            return -1
        elif self.chr_id > other.chr_id:
            return 1
        elif self.start < other.start:
            return -1
        elif self.start > other.start:
            return 1
        elif self.end < other.end:
            return -1
        elif self.end > other.end:
            return 1
        else:
            return 0

    def __lt__(self, other):
        return (self.__cmp__(other) < 0)

    def __le__(self, other):
        return (self.__cmp__(other) <= 0)

    def __gt__(self, other):
        return (self.__cmp__(other) > 0)

    def __ge__(self, other):
        return (self.__cmp__(other) >= 0)

    def __eq__(self, other):
        return (self.__cmp__(other) == 0)

    def __str__(self):
        return f'{self.chr_name}, {self.start}, {self.end}'

    def __repr__(self):
        return self.__str__()


ROOT = os.path.dirname(os.path.abspath(__file__))


filepath = os.path.join(ROOT, 'marker_region_TCGA_6WBC.txt')
print(f'Parsing {filepath}...')
with open(filepath, 'r') as f:
    lines = f.readlines()[1:]
pairs = [
    [{'placenta'}, 'tcga-placenta-sorted.txt'],
    [{'READ', 'STAD', 'BRCA', 'LUSC', 'BLCA', 'LUAD', 'COAD', 'GBM', 'UCEC'}, 'tcga-breast-cancer-sorted.txt'],
]
for tissues, dest_filename in pairs:
    regions = list()
    for line in lines:
        el = line.rstrip().split()
        if (el[0] in tissues) and (int(el[4]) >= 4) and (int(el[5]) >= 150):
            chr_name = el[1]
            if (len(chr_name) > 3) and (chr_name[:3] == 'chr'):
                chr_name = chr_name[3:]
            regions.append(Region(chr_name, el[2], el[3]))
    regions = list(sorted(regions))
    with open(os.path.join(ROOT, dest_filename), 'w') as f:
        for region in regions:
            if not region.chr_id in {23, 24, 25}:
                f.write(f'{region}\n')


for filename in ['marker_region_TCGA_6WBC.txt', 'marker_region_18Normals.txt']:
    filepath = os.path.join(ROOT, filename)
    if not filepath.endswith('txt'):
        continue
    if filepath.endswith('-sorted.txt'):
        continue
    print(f'Parsing {filepath}...')
    with open(filepath, 'r') as f:
        lines = f.readlines()[1:]
    regions = list()
    for line in lines:
        el = line.rstrip().split()
        regions.append(Region(el[1], el[2], el[3]))
    regions = list(sorted(regions))
    new_filepath = os.path.splitext(filepath)[0] + '-sorted.txt'
    with open(new_filepath, 'w') as f:
        for region in regions:
            if not region.chr_id in {23, 24, 25}:
                f.write(f'{region}\n')
