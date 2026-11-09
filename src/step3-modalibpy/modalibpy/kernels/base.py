from typing import Optional

import numpy as np
import scipy.optimize

from sklearn.base import BaseEstimator, TransformerMixin


class Kernel(BaseEstimator, TransformerMixin):

    @staticmethod
    def safe_log(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
        return np.log(np.clip(x, eps, None))

    @staticmethod
    def safe_exp(x: np.ndarray, eps: float = 1e-8, M: float = 100) -> np.ndarray:
        y = np.copy(x)
        is_neginf = (x < -M)
        y[is_neginf] = 0
        is_zero = np.logical_and(-eps < x, x < eps)
        y[is_zero] = 1
        is_posinf = (x > M)
        y[is_posinf] = np.inf
        valid = np.logical_and(-M < x, x < M)
        valid = np.logical_and(valid, ~is_zero)
        y[valid] = np.exp(x[valid])
        return np.nan_to_num(y)

    @staticmethod
    def infer_bandwidth(D2: np.ndarray) -> float:
        scale = 1. / np.mean(D2)

        idx = np.triu_indices(D2.shape[0], k=1)
        d2 = D2[idx[0], idx[1]]

        def f(x):
            gamma = 10 ** x
            k = Kernel.safe_exp(-gamma * scale * d2)
            return -np.std(k)

        res = scipy.optimize.minimize_scalar(f)
        # res = scipy.optimize.minimize_scalar(f, method='bounded', bounds=(-2, 2))
        gamma = 10 ** res.x
        return gamma * scale
