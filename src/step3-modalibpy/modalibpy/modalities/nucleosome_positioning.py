from typing import Optional

import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData
from ..decorrelate import gc_correction
from ..preprocessing import DirichletMultinomialTransformer
from ..pipeline.fs import FeatureSelector


class NucleosomePositioningModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='chi2', **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.feature_weights: Optional[np.ndarray] = None
        self.mask: Optional[np.ndarray] = None

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('Nucleosome positioning (preprocessed)' in data) and self.cache:
            return None

        X = np.clip(data['Nucleosome score histogram'], 0.2, None)
        X = X / np.sum(X, axis=1)[:, np.newaxis]

        self.save(data, 'Nucleosome positioning (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = self.load(data, 'Nucleosome positioning (preprocessed)')

        return self.load(data, 'Nucleosome positioning (preprocessed)')
