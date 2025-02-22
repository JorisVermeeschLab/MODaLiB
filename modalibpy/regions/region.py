class Region:

    def __init__(self, chr_name: str, start: int, end: int) -> None:
        self.chr_id: int = Region.parse_chr_name(chr_name)
        self.start: int = int(start)
        self.end: int = int(end)

    def overlap(self, other: 'Region') -> float:
        fraction = 0.0
        if other.chr_id == self.chr_id:
            if (self.start <= other.start) and (other.start <= self.end) and (self.end < other.end):
                fraction = float(self.end + 1 - other.start) / float(other.end + 1 - other.start)
            elif (other.start < self.start) and (self.start <= other.end) and (other.end <= self.end):
                fraction = float(other.end + 1 - self.start) / float(other.end + 1 - other.start)
            elif (self.start <= other.start) and (other.end <= self.end):
                fraction = 1.0
            elif (other.start < self.start) and (self.end < other.end):
                fraction = float(self.end + 1 - self.start) / float(other.end + 1 - other.start)
        return fraction

    @staticmethod
    def parse_chr_name(chr_name: str) -> int:
        chr_name = chr_name.upper()
        if chr_name.startswith('CHR'):
            chr_name = chr_name[3:]
        chr_name = chr_name.strip()
        if chr_name == 'X':
            chr_id = 23
        elif chr_name == 'Y':
            chr_id = 24
        elif chr_name in {'M', 'MT'}:
            chr_id = 25
        elif chr_name == 'EBV':
            chr_id = 26
        else:
            chr_id = int(chr_name)
        return chr_id
