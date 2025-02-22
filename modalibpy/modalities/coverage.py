from typing import Optional

import tqdm
import numpy as np
import scipy.signal
from sklearn.preprocessing import RobustScaler, StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold

from .base import BaseModality, MultimodalData
from ..decorrelate import gc_correction
from ..deconvolution.optimizer import Optimizer
from ..preprocessing import DirichletMultinomialTransformer
from ..pipeline.fs import FeatureSelector


class CoverageModality(BaseModality):

    def __init__(self, *args, caff: bool = False, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='rbf', **kwargs)
        self.fs: FeatureSelector = FeatureSelector(n_min=50, n_max=500)
        self.dm_transformer: DirichletMultinomialTransformer = DirichletMultinomialTransformer(gamma=10, max_n_iter=100)
        self.use_caff: bool = caff

    @property
    def key(self) -> str:
        return 'CAFF (preprocessed)' if self.use_caff else 'Coverage profiles (preprocessed)'

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if (self.key in data) and self.cache:
            return None

        X = data['CAFF'] if self.use_caff else data['Coverage']

        # Smoothing using a moving average (window size = 20)
        winsize = 20
        mid_winsize = winsize // 2
        Y = np.cumsum(X, axis=1)
        X[:, mid_winsize:-mid_winsize] = (Y[:, winsize:] - Y[:, :-winsize]) / winsize

        # Correct for the sequencing depth using the median
        medians = np.median(X, axis=1)[:, np.newaxis]
        medians[medians == 0] = 1
        X = X / medians

        self.save(data, self.key, X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = self.load(data, self.key)

        # Feature selection
        #if fit:
        #    X = self.fs.fit_transform(X, y)
        #else:
        #    X = self.fs.transform(X)

        return X


class CAFFModality(CoverageModality):

    def __init__(self, *args, **kwargs) -> None:
        CoverageModality.__init__(self, *args, caff=True, **kwargs)
