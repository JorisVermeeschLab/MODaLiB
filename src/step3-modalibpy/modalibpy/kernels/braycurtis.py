from typing import Optional

import numpy as np
import scipy.spatial

from sklearn.base import BaseEstimator, TransformerMixin

from .base import Kernel


class BrayCurtisKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'RBFKernel':
        self.X_train = X
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        D = scipy.spatial.distance.cdist(X, self.X_train, metric='braycurtis')
        return 1 - D
