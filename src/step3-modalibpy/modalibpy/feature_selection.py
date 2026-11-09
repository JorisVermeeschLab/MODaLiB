from typing import Optional

import tqdm
import numpy as np
import scipy.stats
from sklearn.linear_model import LogisticRegression
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import KBinsDiscretizer
from pyitlib import discrete_random_variable as drv


class CMISelector(BaseEstimator, TransformerMixin):

    def __init__(self, k: int):
        self.k: int = k
        self.mask: Optional[np.ndarray] = None

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None, c: Optional[np.ndarray] = None) -> 'CMISelector':

        if c is None:
            cmi = mutual_info_classif(X, y)
        else:
            disc = KBinsDiscretizer(n_bins=5, encode='ordinal', strategy='uniform')
            X = disc.fit_transform(X)
            cmi = drv.information_mutual_conditional(X.T.astype(int), y.astype(int), c.astype(int), cartesian_product=True)

        idx = np.argsort(cmi)[-self.k:]
        self.mask = np.zeros(X.shape[1], dtype=bool)
        self.mask[idx] = True

        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        if self.mask is not None:
            X = X[:, self.mask]
        return X


class CausalFeatureSelector(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.mask: Optional[np.ndarray] = None

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'CausalFeatureSelector':

        residuals_y = np.zeros(X.shape, dtype=int)
        residuals_x = np.zeros_like(X)
        for j in tqdm.tqdm(range(X.shape[1])):

            model = LogisticRegression()
            model.fit(X[:, j, np.newaxis], y)
            y_hat = model.predict(X[:, j, np.newaxis])
            residuals_y[:, j] = (y_hat != y).astype(int)

            for c in np.unique(y):
                mask = (y == c)
                residuals_x[mask, j] = (np.mean(X[mask, j]) - X[mask, j])

        mis_x = np.squeeze([mutual_info_classif(X[:, j, np.newaxis], residuals_y[:, j]) for j in range(X.shape[1])])
        mis_y = mutual_info_classif(residuals_x, y)

        import matplotlib.pyplot as plt
        plt.scatter(mis_x, mis_y, alpha=0.3)
        plt.show()

        self.mask = (mis_x < mis_y)

        print(np.mean(self.mask))

        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        if self.mask is not None:
            X = X[:, self.mask]
        return X
