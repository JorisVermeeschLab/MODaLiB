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
        self.dm_transformer: DirichletMultinomialTransformer = DirichletMultinomialTransformer(gamma=10, max_n_iter=2000)
        self.fs: FeatureSelector = FeatureSelector(n_min=50)
        self.mask: Optional[np.ndarray] = None

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('Nucleosome positioning (preprocessed)' in data) and self.cache:
            return None

        X1 = data['Nucleosome positioning']
        n = data['Nucleosome score counts']

        # Smoothing using a moving average (window size = 20)
        #winsize = 20
        #mid_winsize = winsize // 2
        #Y = np.cumsum(X1, axis=1)
        #X1[:, mid_winsize:-mid_winsize] = (Y[:, winsize:] - Y[:, :-winsize]) / winsize
        #X1 = np.clip(X1, 0.01, 0.99)
        #X, n = X1[:, ::winsize], X1[:, ::winsize]

        # Correcting with baseline positioning at CTCF-binding sites
        #X2 = data['Nucleosome score histogram (CTCF)']
        #X2 = X2 / np.sum(X2, axis=1)[:, np.newaxis]
        #x2 = np.sum(np.arange(X2.shape[1])[np.newaxis, :] * X2, axis=1)[:, np.newaxis]
        #X = np.log(X1 / x2)

        X = np.clip(data['Nucleosome score histogram'], 0.2, None)

        if self.probabilistic:
            if fit:
                X = self.dm_transformer.fit_transform(X)
            else:
                X = self.dm_transformer.transform(X)
        else:
            X = X / np.sum(X, axis=1)[:, np.newaxis]

        # Remove low-coverage bins
        #n = data['Nucleosome score counts']
        #if fit:
        #    self.mask = (np.mean(n, axis=0) > 200)
        #X = X[:, self.mask]

        self.save(data, 'Nucleosome positioning (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = self.load(data, 'Nucleosome positioning (preprocessed)')

        # Feature selection
        #if fit:
        #    X = self.fs.fit_transform(X, y)
        #else:
        #    X = self.fs.transform(X)

        return self.load(data, 'Nucleosome positioning (preprocessed)')
