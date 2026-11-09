from typing import Optional

import numpy as np
import scipy.spatial

from sklearn.base import BaseEstimator, TransformerMixin

from .base import Kernel


class RBFKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'RBFKernel':
        self.X_train = X
        D2 = scipy.spatial.distance.cdist(X, X, metric='sqeuclidean')
        #self.gamma = Kernel.infer_bandwidth(D2)
        self.gamma = 1. / (X.shape[1] * X.var())
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        D2 = scipy.spatial.distance.cdist(X, self.X_train, metric='sqeuclidean')
        return Kernel.safe_exp(-self.gamma * D2)
