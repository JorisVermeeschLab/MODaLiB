from dataclasses import dataclass


@dataclass(init=True, eq=True)
class Fragment:
    start_chr_id: int = 0
    start_pos: int = 0
    end_chr_id: int = 0
    end_pos: int = 0
    size: int = 0
    n_methylated: int = 0
    n_unmethylated: int = 0
    snp_count: int = 0
    gc_count: int = 0
    s_count: int = 0
    n_inserted: int = 0
    n_deleted: int = 0


"""
class Fragment:

    def __init__(
        self,
        start_chr_id: int = 0,
        start_pos: int = 0,
        end_chr_id: int = 0,
        end_pos: int = 0,
        size: int = 0,
        n_methylated: int = 0,
        n_unmethylated: int = 0,
        snp_count: int = 0,
        gc_count: int = 0,
        s_count: int = 0,
        n_inserted: int = 0,
        n_deleted: int = 0
    ):
        self.start_chr_id: int = start_chr_id
        self.start_pos: int = start_pos
        self.end_chr_id: int = end_chr_id
        self.end_pos: int = end_pos
        self.size: int = size
        self.n_methylated: int = n_methylated
        self.n_unmethylated: int = n_unmethylated
        self.snp_count: int = snp_count
        self.gc_count: int = gc_count
        self.s_count: int = s_count
        self.n_inserted: int = n_inserted
        self.n_deleted: int = n_deleted

    def __repr__(self) -> str:
        s = f'Fragment(start={self.start_chr_id}-{self.start_pos}, end={self.end_chr_id}-{self.end_pos}, size={self.size}, '
        s += f'meth={self.n_methylated}/{self.n_unmethylated}, snp={self.snp_count})'
        return s

    def __str__(self) -> str:
        return self.__repr__()
"""