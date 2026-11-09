from typing import Optional

import numpy as np
import scipy.spatial

from sklearn.base import BaseEstimator, TransformerMixin

from .base import Kernel


euclidean_dist_matrix = scipy.spatial.distance.cdist


class ANOVAKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'RBFKernel':
        self.X_train = X
        #D2 = scipy.spatial.distance.cdist(X, X, metric='sqeuclidean')
        #self.gamma = Kernel.infer_bandwidth(D2)
        #return self
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:

        data_1 = X
        data_2 = self.X_train


        self._q = 0.1
        kernel = np.ones((data_1.shape[0], data_2.shape[0]))

        for d in range(data_1.shape[1]):
            column_1 = data_1[:, d].reshape(-1, 1)
            column_2 = data_2[:, d].reshape(-1, 1)
            kernel *= (1-self._q ** 2) / \
                      (2.*(1. - 2.*self._q *np.cos(column_1 - column_2.T) + self._q ** 2))

        return kernel
