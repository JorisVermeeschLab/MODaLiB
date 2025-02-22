from typing import Optional

import numpy as np

from .base import BaseModality, MultimodalData
from ..preprocessing import DirichletMultinomialTransformer


class MutationalSignatureModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='chi2', **kwargs)
        self.dm_transformer: DirichletMultinomialTransformer = DirichletMultinomialTransformer()

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if 'Mutational signature (preprocessed)' in data:
            return None

        counts = data['Mutations']
        counts = np.clip(counts, 0.1, None)

        if self.probabilistic:
            Y = np.log1p(np.repeat(data['Library size'][:, np.newaxis], counts.shape[1], axis=1))
            if fit:
                X = self.dm_transformer.fit_transform(counts, y=Y)
            else:
                X = self.dm_transformer.transform(counts, y=Y)
        else:
            X = counts / np.sum(counts, axis=1)[:, np.newaxis]
        self.save(data, 'Mutational signature (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        return self.load(data, 'Mutational signature (preprocessed)')
