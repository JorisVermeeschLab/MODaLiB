from typing import List, Dict

import tqdm
import numpy as np
import scipy.stats
import seaborn
import matplotlib.pyplot as plt


class CellRanking:

    def __init__(
            self,
            X: np.ndarray,
            labels: List[str],
            R: np.ndarray,
            cell_types: List[str],
            tissues: List[str],
            compartments: List[str]
    ):
        self.X: np.ndarray = X
        self.labels: np.ndarray = np.asarray(labels, dtype=object)
        self.R: np.ndarray = R
        self.cell_types: np.ndarray = np.asarray(cell_types, dtype=object)
        self.tissues: np.ndarray = np.asarray(tissues, dtype=object)
        self.compartments: np.ndarray = np.asarray(compartments, dtype=object)

        n = len(X)
        m = len(R)
        p = X.shape[1]
        assert len(labels) == n
        assert len(cell_types) == m

        self.C = np.zeros((n, m))
        for i in tqdm.tqdm(range(len(X)), desc='Cell ranking'):
            for j in range(len(R)):
                self.C[i, j] = scipy.stats.pearsonr(X[i, :], R[j, :]).statistic

    def plot(self, palette: Dict[str, str]) -> None:
        n = len(self.X)
        m = len(self.R)
        ranks = m - np.argsort(np.abs(self.C), axis=1)
        #ranks = np.abs(self.C)

        unique_compartments = ['immune', 'endothelial', 'stromal', 'epithelial']
        ranks_per_compartment = []
        for compartment in unique_compartments:
            mask = (self.compartments == compartment)
            ranks_per_compartment.append(np.mean(np.abs(ranks[:, mask]), axis=1))
        ranks_per_compartment = np.asarray(ranks_per_compartment).T

        print(ranks_per_compartment.shape)

        seaborn.reset_orig()
        xs, ys, hue = [], [], []
        for i in range(n):
            for j in range(len(unique_compartments)):
                xs.append(ranks_per_compartment[i, j])
                ys.append(unique_compartments[j])
                hue.append(self.labels[i])
        seaborn.violinplot(x=xs, y=ys, hue=hue, palette=palette, saturation=1, fill=False, inner='quart')

