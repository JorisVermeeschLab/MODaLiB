from typing import Optional

import numpy as np
import scipy.spatial

from sklearn.base import BaseEstimator, TransformerMixin


class RBFKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'RBFKernel':
        self.X_train = X

        D2 = scipy.spatial.distance.cdist(X, X, metric='sqeuclidean')

        scale = 1. / np.mean(D2)

        idx = np.triu_indices(D2.shape[0], k=1)
        d2 = D2[idx[0], idx[1]]

        def f(x):
            gamma = 10 ** x
            k = np.exp(-gamma * scale * d2)
            return -np.std(k)

        res = scipy.optimize.minimize_scalar(f)
        # res = scipy.optimize.minimize_scalar(f, method='bounded', bounds=(-2, 2))
        gamma = 10 ** res.x

        self.gamma = gamma * scale
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        D2 = scipy.spatial.distance.cdist(X, self.X_train, metric='sqeuclidean')
        return np.exp(-self.gamma * D2)
