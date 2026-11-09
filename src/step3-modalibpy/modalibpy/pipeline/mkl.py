from typing import Optional, List, Tuple

import numpy as np
import scipy.stats
from sklearn.decomposition import PCA, KernelPCA
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
import torch

from ..modalities.base import MultimodalData
from ..modalities import MODALITY_CLASSES, Modalities
from ..kernels.base import Kernel


class KernelMatrixExtractor:

    def __init__(self, modality_enums: List[Modalities], mod_options: dict=dict()):
        self.modality_enums: List[Modalities] = modality_enums
        self.mod_options = mod_options
        self.modalities = []
        for modality_enum in self.modality_enums:
            modality = MODALITY_CLASSES[modality_enum](**self.mod_options)
            self.modalities.append(modality)

    def preprocess(self, data: MultimodalData, fit: bool=True) -> None:
        kernel_matrices = []
        for modality in self.modalities:
            if fit:
                kernel_matrices.append(modality.preprocess(data, fit=True))
            else:
                kernel_matrices.append(modality.preprocess(data, fit=False))

    def compute_kernel_matrices(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool=True) -> np.ndarray:
        kernel_matrices = []
        for modality in self.modalities:
            if fit:
                kernel_matrices.append(modality.fit_transform(data, y=y))
            else:
                kernel_matrices.append(modality.transform(data))
        return np.asarray(kernel_matrices)

    def fit(self, data: MultimodalData, y: Optional[np.ndarray] = None) -> 'KernelMatrixExtractor':
        return self

    def fit_transform(self, data: MultimodalData, y: Optional[np.ndarray] = None) -> np.ndarray:
        self.fit(data, y=y)
        return self.compute_kernel_matrices(data, y=y, fit=True)

    def transform(self, data: MultimodalData) -> np.ndarray:
        return self.compute_kernel_matrices(data, fit=False)


class MultipleKernelLearning:

    def __init__(self, supervised: bool = False):
        self.supervised: bool = supervised
        self.gamma: float = 1.0
        self.weights: Optional[np.ndarray] = None

    @staticmethod
    def compute_ideal_kernel_matrix(y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        y = y[:, np.newaxis]
        n = len(y)
        K = np.zeros((n, n))
        W = np.zeros((n, n))

        # Both samples are controls -> Should be similar
        #mask = np.logical_and(y == 0, y.T == 0)
        mask = (y == y.T)
        K[mask] = 1
        W[mask] = 1. / int(np.sum(mask))

        # One sample is a control and the other one is a cancer case -> Should be dissimilar
        mask = (y != y.T)
        K[mask] = 0
        W[mask] = 1. / int(np.sum(mask))

        # Both samples are cancer cases -> Inconclusive
        #mask = np.logical_and(y == 1, y.T == 1)
        #K[mask] = 1
        #W[mask] = 0  # 1. / int(np.sum(mask))  #0  # These sample pairs will be ignored

        W /= np.sum(W)
        return K, W

    @staticmethod
    def infer_weights_unsupervised(kernel_matrices: np.ndarray) -> np.ndarray:

        X = []
        for k, K in enumerate(kernel_matrices):
            X.append(KernelPCA(n_components=len(K)).fit_transform(K))
        X = np.concatenate(X, axis=1)

        C = np.abs(np.cov(X.T))
        weights = np.mean(C, axis=0)

        weights = weights.reshape(len(kernel_matrices), -1)
        weights = np.mean(weights, axis=1)
        weights += 1e-10
        weights /= np.sum(weights)
        return weights

    @staticmethod
    def weighted_pearson_corrcoef(x, y, w):
        x = np.argsort(x)
        y = np.argsort(y)
        mean_x = np.sum(x * w) / np.sum(w)
        mean_y = np.sum(y * w) / np.sum(w)
        var_x = np.sum(w * (x - mean_x) * (x - mean_x)) / np.sum(w)
        var_y = np.sum(w * (y - mean_y) * (y - mean_y)) / np.sum(w)
        cov = np.sum(w * (x - mean_x) * (y - mean_y)) / np.sum(w)
        return cov / np.sqrt(var_x * var_y)

    @staticmethod
    def conditional_correlation(x: np.ndarray, y: np.ndarray, z: np.ndarray, weights: np.ndarray) -> float:

        # Regress x on z
        model = IsotonicRegression(increasing=True, y_min=np.min(x), y_max=np.max(x))
        model.fit(z, x, sample_weight=weights)
        x_hat = model.predict(z)
        mask = np.isnan(x_hat)
        if np.any(mask):
            model = IsotonicRegression(increasing=True, y_min=np.min(x), y_max=np.max(x))
            model.fit(z, x)
            x_hat[mask] = model.predict(z[mask])
        x_res = x - x_hat
        assert not np.any(np.isnan(x))
        assert not np.any(np.isnan(z))
        assert not np.any(np.isnan(x_hat))

        # Regress y on z
        model = LogisticRegression()
        model.fit(z[:, np.newaxis], y, sample_weight=weights)
        y_hat = model.predict_proba(z[:, np.newaxis])[:, 1]
        y_res = y - y_hat
        assert not np.any(np.isnan(y_res))

        #return scipy.stats.pearsonr(x_res, y_res).statistic
        return MultipleKernelLearning.weighted_pearson_corrcoef(x_res, y_res, weights)

    @staticmethod
    def infer_weights(Ks: np.ndarray, y: np.ndarray) -> np.ndarray:

        if len(Ks) == 1:
            return np.asarray([1.])

        n_modalities = len(Ks)
        K_y, W_y = MultipleKernelLearning.compute_ideal_kernel_matrix(y)

        weights = []
        for p in range(len(Ks)):
            mask = np.ones(len(Ks), dtype=bool)
            mask[p] = False
            #K_mean = np.mean(Ks[mask, :, :], axis=0)
            K_mean = Kernel.safe_exp(np.mean(Kernel.safe_log(Ks[mask, :, :]), axis=0))
            weights.append(MultipleKernelLearning.conditional_correlation(
                Ks[p, :, :].flatten(),
                K_y.flatten(),
                K_mean.flatten(),
                W_y.flatten()
            ))
        weights = np.asarray(weights)
        weights = np.clip(weights, 1e-10, None)
        weights /= np.sum(weights)

        return weights

    @staticmethod
    def preprocess(Ks: np.ndarray) -> np.ndarray:
        #n_modalities, n, m = Ks.shape
        #return np.sqrt(np.einsum('ikl,jkl->ijkl', Ks, Ks)).reshape(n_modalities * n_modalities, n, m)
        return Ks

    def fit(self, Ks: np.ndarray, y: np.ndarray) -> 'MultipleKernelLearning':
        Ks = MultipleKernelLearning.preprocess(Ks)
        if self.supervised:
            self.weights = MultipleKernelLearning.infer_weights(Ks, y)
        else:
            self.weights = MultipleKernelLearning.infer_weights_unsupervised(Ks)
        return self

    def fit_transform(self, Ks: np.ndarray, y: np.ndarray) -> np.ndarray:
        self.fit(Ks, y)
        return self.transform(Ks)

    def transform(self, Ks: np.ndarray) -> np.ndarray:
        Ks = MultipleKernelLearning.preprocess(Ks)
        Ks = Kernel.safe_log(Ks)
        Ks = Ks * self.weights[:, np.newaxis, np.newaxis]
        K = Ks.sum(axis=0)
        K = Kernel.safe_exp(K)
        return K
