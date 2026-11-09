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
        X = counts / np.sum(counts, axis=1)[:, np.newaxis]

        self.save(data, 'Fragment length counts (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        return self.load(data, 'Fragment length counts (preprocessed)')
