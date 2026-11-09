from typing import List, Generator

import pandas as pd


class Region:

    def __init__(self, chr_name: str, start: int, end: int, name: str):
        self.chr_name: str = chr_name
        self.start: int = start
        self.end: int = end
        self.name: str = name

    def __contains__(self, other: 'Region') -> bool:
        if self.chr_name != other.chr_name:
            return False
        return (self.start <= other.start) and (other.end <= self.end)

    def __str__(self) -> str:
        return f'{self.chr_name}:{self.start}-{self.end} ({self.name})'

    def __repr__(self) -> str:
        return self.__str__()


class Regions:

    def __init__(self):
        self.regions: List[Region] = []

    def add(self, region: Region) -> None:
        self.regions.append(region)

    def __getitem__(self, i: int) -> Region:
        return self.regions[i]

    def __iter__(self) -> Generator[None, None, Region]:
        for region in self.regions:
            yield region

    def __len__(self) -> int:
        return len(self.regions)

    @staticmethod
    def load(filepath: str) -> 'Regions':
        regions = Regions()
        with open(filepath, 'r') as f:
            for line in f:
                elements = line.rstrip().split(',')
                if len(elements) >= 4:
                    chr_id = int(elements[0])
                    if chr_id == 23:
                        chr_id = 'X'
                    elif chr_id == 24:
                        chr_id = 'chrY'
                    elif chr_id == 25:
                        chr_id = 'MT'
                    elif chr_id == 26:
                        chr_id = 'EBV'
                    else:
                        chr_id = str(chr_id)
                    regions.add(Region(chr_id, int(elements[1]), int(elements[2]), elements[3]))
        return regions

    @staticmethod
    def load_cytobands(filepath: str) -> 'Regions':
        df = pd.read_csv(filepath, sep='\t', header=None, names=['chr', 'start', 'end', 'band', 'stain'])
        df = df[~df['band'].isna()]
        regions = Regions()
        for c, start, end, band in zip(df['chr'], df['start'], df['end'], df['band']):
            chr_name = c.replace('chr', '')
            regions.add(Region(chr_name, int(start), int(end), chr_name + band))
        return regions
