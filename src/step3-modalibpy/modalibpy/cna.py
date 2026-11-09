from typing import Dict, Set

import numpy as np
from sklearn.preprocessing import RobustScaler

from modalibpy.region import Region, Regions


class CoverageData(object):

    def __init__(self, X: np.ndarray, bins: Regions, cytobands: Regions, genes: Regions):
        assert len(X.shape) == 2
        self.X: np.ndarray = X
        self.bins: Regions = bins
        self.cytobands: Regions = cytobands
        self.cytoband_dict: Dict[str, Region] = {cytoband.name: cytoband for cytoband in cytobands}
        self.genes: Regions = genes

    def __getitem__(self, cytoband_name: str) -> np.ndarray:
        idx = self.get_chromosome_region(cytoband_name)
        if len(idx) == 0:
            x = np.zeros(self.n_samples)
        else:
            x = np.sum(self.X[:, idx], axis=1)
        return x

    def _get_chromosome_band(self, cytoband: Region) -> Set[int]:
        idx = set()
        for i, r in enumerate(self.bins):
            if r in cytoband:
                idx.add(i)
        return idx

    def get_chromosome_region(self, cytoband_name: str) -> np.ndarray:
        if '.' in cytoband_name:
            cytoband = self.cytoband_dict[cytoband_name]
            idx = self._get_chromosome_band(cytoband)
        else:
            idx = set()
            for key in self.cytoband_dict:
                if key.split('.')[0].startswith(cytoband_name):
                    idx = idx.union(self._get_chromosome_band(self.cytoband_dict[key]))
        idx = np.asarray(list(idx))
        return np.sort(idx)

    def get_gene_set_(self, cytoband: Region) -> Set[str]:
        gene_names = set()
        for gene in self.genes:
            if gene in cytoband:
                gene_names.add(gene.name)
        return gene_names

    def get_gene_set(self, cytoband_name: str) -> Set[str]:
        gene_names = set()
        if '.' in cytoband_name:
            cytoband = self.cytoband_dict[cytoband_name]
            gene_names.update(self.get_gene_set_(cytoband))
        else:
            for key in self.cytoband_dict:
                if key.split('.')[0].startswith(cytoband_name):
                    gene_names.update(self.get_gene_set_(self.cytoband_dict[key]))
        return gene_names

    @property
    def n_samples(self):
        return len(self.X)
