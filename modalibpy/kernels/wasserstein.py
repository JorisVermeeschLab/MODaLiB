from typing import Optional

import numpy as np
import scipy.spatial
import scipy.special

from sklearn.base import BaseEstimator, TransformerMixin

from .base import Kernel


class WassersteinKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    @staticmethod
    def pairwise_distances(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        #w = np.linspace(0, 1, X.shape[1] + 3)[1:-1]
        ##w = scipy.special.logit(w)
        #w = np.abs(w[1:] - w[:-1])
        #w /= np.sum(w)
        #return scipy.spatial.distance.cdist(X, Y, metric='cityblock', w=w) * X.shape[1]
        return scipy.spatial.distance.cdist(X, Y, metric='cityblock')

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'WassersteinKernel':
        X = X / np.sum(X, axis=1)[:, np.newaxis]
        self.X_train = X
        X = np.cumsum(X, axis=1)
        D2 = WassersteinKernel.pairwise_distances(X, X)
        self.gamma = Kernel.infer_bandwidth(D2)
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        X = X / np.sum(X, axis=1)[:, np.newaxis]
        D2 = WassersteinKernel.pairwise_distances(
            np.cumsum(X, axis=1),
            np.cumsum(self.X_train, axis=1)
        )
        return Kernel.safe_exp(-self.gamma * D2)


class GaussianWassersteinKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.X_train: Optional[np.ndarray] = None
        self.gamma: float = 1

    @staticmethod
    def cdist(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        X1 = X[:, :, 0]
        X2 = X[:, :, 1]
        Y1 = Y[:, :, 0]
        Y2 = Y[:, :, 1]
        mu1 = X1
        mu2 = Y1
        var1 = X2 - np.square(X1)
        var2 = Y2 - np.square(Y1)
        D = scipy.spatial.distance.cdist(mu1, mu2, metric='sqeuclidean')
        D += np.sum(var1, axis=1)[:, np.newaxis]
        D += np.sum(var2, axis=1)[np.newaxis, :]
        D -= 2 * np.dot(np.sqrt(var1), np.sqrt(var2.T))
        return D

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'WassersteinKernel':
        self.X_train = X
        D2 = GaussianWassersteinKernel.cdist(X, X)
        self.gamma = Kernel.infer_bandwidth(D2)
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        D2 = GaussianWassersteinKernel.cdist(X, self.X_train)
        return Kernel.safe_exp(-self.gamma * D2)
