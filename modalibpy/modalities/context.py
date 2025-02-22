from typing import Optional

import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import DirichletMultinomialTransformer


class TrinucleotideContext(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='chi2', **kwargs)
        self.dm_transformer: DirichletMultinomialTransformer = DirichletMultinomialTransformer()
        self.scaler: RobustScaler = RobustScaler()

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if 'Trinucleotide context (preprocessed)' in data:
            return None

        counts = data['Trinucleotide frequencies']
        counts = np.clip(counts, 0.1, None)

        if self.probabilistic:
            Y = np.log1p(np.repeat(data['Library size'][:, np.newaxis], counts.shape[1], axis=1))
            if fit:
                X = self.dm_transformer.fit_transform(counts, y=Y)
            else:
                X = self.dm_transformer.transform(counts, y=Y)
        else:
            X = counts / np.sum(counts, axis=1)[:, np.newaxis]

        self.save(data, 'Trinucleotide context (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        return self.load(data, 'Trinucleotide context (preprocessed)')
