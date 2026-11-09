from typing import Optional

import numpy as np
import scipy.stats
import sklearn.feature_selection


class FeatureSelector:

    def __init__(self, n_min: int = 25, n_max: int = np.inf):
        self.n_min: int = n_min
        self.n_max: int = n_max
        self.selected_idx: Optional[np.ndarray] = None

    def fit(self, X: np.ndarray, y: np.ndarray = None) -> 'FeatureSelector':
        if y is not None:
            _, p_values = sklearn.feature_selection.f_classif(X, y)
            #p_values = np.asarray([scipy.stats.ttest_ind(X[y == 0, k], X[y == 1, k]).pvalue for k in range(X.shape[1])], dtype=float)
            #p_values = np.asarray([scipy.stats.ks_2samp(X[y == 0, k], X[y == 1, k]).pvalue for k in range(X.shape[1])], dtype=float)
            sig_cutoff = 0.001 / len(p_values)
            mask = (p_values < sig_cutoff)
            idx = np.argsort(p_values)
            n_max = len(idx) if (self.n_max > len(idx)) else self.n_max
            mask[idx[n_max:]] = False
            if np.sum(mask) < self.n_min:
                self.selected_idx = idx[:self.n_min]
            else:
                self.selected_idx = np.where(mask)[0]
        else:
            self.selected_idx = np.arange(X.shape[1])
        return self

    def fit_transform(self, X: np.ndarray, y: np.ndarray = None):
        self.fit(X, y)
        return self.transform(X)

    def transform(self, X: np.ndarray) -> np.ndarray:
        return X[:, self.selected_idx]
