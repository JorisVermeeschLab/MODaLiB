from typing import Optional

import numpy as np
import scipy.stats
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA

from .base import BaseModality, MultimodalData
from ..end_motifs import cgnngc_ratios
from ..preprocessing import BetaBinomialTransformer
from ..decorrelate import gc_correction
from ..pipeline.fs import FeatureSelector


class CpGMethylationModality(BaseModality):

    def __init__(self, *args, fmr: bool = False, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='rbf', **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.bb_transformer: BetaBinomialTransformer = BetaBinomialTransformer(gamma=0.01, lr=1e-7, max_n_iter=50)
        self.feature_weights: Optional[np.ndarray] = None
        self.fs: FeatureSelector = FeatureSelector(n_min=500)
        self.use_fmr: bool = fmr

    @property
    def key(self) -> str:
        return 'FMR (preprocessed)' if self.use_fmr else 'CpG methylation (preprocessed)'

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if (self.key in data) and self.cache:
            return None

        X = data['Fully methylated fragment ratio'] if self.use_fmr else data['CpG Methylation (first moment)']
        n = np.clip(data['Coverage'] if self.use_fmr else data['CpG count'], 0.01, None) * 10
        
        y_cgn_ncg = cgnngc_ratios(data['End-motif probability (5 prime)'])
        Y = np.repeat(y_cgn_ncg[:, np.newaxis], X.shape[1], axis=1)

        if self.probabilistic:
            if fit:
                X = self.bb_transformer.fit_transform(n, X * n, y=Y)
            else:
                X = self.bb_transformer.transform(n, X * n, y=Y)

        self.save(data, self.key, X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = self.load(data, self.key)
        n = np.clip(data['Coverage'] if self.use_fmr else data['CpG count'], 0.01, None)

        # Feature selection
        if fit:
            X = self.fs.fit_transform(X, y)
        else:
            X = self.fs.transform(X)

        return X


class FullyMethylatedFragModality(CpGMethylationModality):

    def __init__(self, *args, **kwargs) -> None:
        CpGMethylationModality.__init__(self, *args, fmr=True, **kwargs)
