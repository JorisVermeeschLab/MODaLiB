from typing import Optional

import numpy as np
import scipy.spatial

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics.pairwise import additive_chi2_kernel

from .base import Kernel


class ChiSquaredKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'RBFKernel':
        self.X_train = X
        D2 = -additive_chi2_kernel(X, X)
        self.gamma = Kernel.infer_bandwidth(D2)
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        D2 = -additive_chi2_kernel(X, self.X_train)
        return Kernel.safe_exp(-self.gamma * D2)
