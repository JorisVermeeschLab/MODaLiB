from typing import Optional

import numpy as np
import scipy.stats
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import DirichletMultinomialTransformer
from ..end_motifs import *
from .jaggedness import JaggednessModality


def average_frag_lengths(data: MultimodalData) -> np.ndarray:
    W = data['Fragment length counts']
    W = W / np.sum(W, axis=1)[:, np.newaxis]
    lengths = np.arange(W.shape[1])
    return np.sum(lengths[np.newaxis, :] * W, axis=1)


class EndMotifFrequencyModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='chi2', **kwargs)
        self.dm_transformer: DirichletMultinomialTransformer = DirichletMultinomialTransformer(gamma=10, max_n_iter=1000, lr=0.01)
        self.scaler: RobustScaler = RobustScaler()

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('End motif frequencies (preprocessed)' in data) and self.cache:
            return None

        counts = data['End-motif probability (5 prime)']
        counts = np.clip(counts, 0.1, None) 

        """
        if self.probabilistic:  # TODO
            y = average_frag_lengths(data)
            #y = np.mean(data['GC-content'], axis=1)
            #y = JaggednessModality.compute_jaggedness(data)
            #y = semfar_scores(counts)
            Y = np.repeat(y[:, np.newaxis], counts.shape[1], axis=1)
            Y = RobustScaler().fit_transform(Y)
            if fit:
                X = self.dm_transformer.fit_transform(counts, y=Y)
            else:
                X = self.dm_transformer.transform(counts, y=Y)
        else:
            X = counts / np.sum(counts, axis=1)[:, np.newaxis]
        """

        X = counts / np.sum(counts, axis=1)[:, np.newaxis]

        self.save(data, 'End motif frequencies (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        X = self.load(data, 'End motif frequencies (preprocessed)')
        return X
