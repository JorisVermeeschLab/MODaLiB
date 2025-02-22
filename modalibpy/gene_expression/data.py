from typing import List, Tuple

import numpy as np


class GeneExpressionData(object):

    def __init__(self, X: np.ndarray, gene_names: List[str]):
        self.X: np.ndarray = X
        self.gene_names: List[str] = gene_names

    def reorder(self, gene_names: List[str]) -> 'GeneExpressionData':
        gene_dict = {gene_name: i for i, gene_name in enumerate(self.gene_names)}
        idx = np.asarray([gene_dict[gene_name] for gene_name in gene_names], dtype=int)
        return GeneExpressionData(self.X[:, idx], gene_names)

    def match(self, other: 'GeneExpressionData') -> Tuple['GeneExpressionData', 'GeneExpressionData']:
        gene_names = list(set(self.gene_names).intersection(other.gene_names))
        return self.reorder(gene_names), other.reorder(gene_names)

