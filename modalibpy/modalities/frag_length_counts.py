from typing import Optional

import numpy as np
from sklearn.preprocessing import StandardScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import DirichletMultinomialTransformer
from ..end_motifs import *
from .jaggedness import JaggednessModality


class FragLengthCountModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='rbf', **kwargs)
        self.dm_transformer: DirichletMultinomialTransformer = DirichletMultinomialTransformer(gamma=10, max_n_iter=1000, lr=0.005)

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('Fragment length counts (preprocessed)' in data) and self.cache:
            return None

        counts = data['Fragment length counts']
        counts = counts[:, 40:] * 10

        """
        if self.probabilistic:  # TODO
            #y = semfar_scores(data['End-motif probability (5 prime)'])
            #y = JaggednessModality.compute_jaggedness(data)
            #y = np.mean(data['Mappability'], axis=1)
            y = semfar_scores(data['End-motif probability (5 prime)'])
            #y = np.mean(data['GC-content'], axis=1)
            #y = data['Library size']
            Y = np.repeat(y[:, np.newaxis], counts.shape[1], axis=1)
            Y = StandardScaler().fit_transform(Y)
            #Y = None
            if fit:
                X = self.dm_transformer.fit_transform(counts, y=Y)
            else:
                X = self.dm_transformer.transform(counts, y=Y)
        else:
            X = counts / np.sum(counts, axis=1)[:, np.newaxis]
        """
        X = counts / np.sum(counts, axis=1)[:, np.newaxis]

        self.save(data, 'Fragment length counts (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        return self.load(data, 'Fragment length counts (preprocessed)')
