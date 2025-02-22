import copy
import enum
from typing import Tuple, Optional, List

import tqdm
import numpy as np
import scipy.spatial
import scipy.stats
import scipy.optimize
import scipy.special
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import pairwise_kernels, chi2_kernel, additive_chi2_kernel
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer
from sklearn.feature_selection import f_classif, SelectKBest
from sklearn.decomposition import KernelPCA

from modalibpy.kernels import RBFKernel
from modalibpy.distributions import BinomialDistribution, BetaBinomialModel
from .decorrelate import Decorrelate
from ..modalities.base import MultimodalData
from ..modalities import MODALITY_CLASSES, Modalities
from .scaler import Scaler
from .reweight import Reweight


class FeaturesExtractor:

    def __init__(self, modality_enums: List[Modalities], mod_options: dict=dict()):
        self.modality_enums: List[Modalities] = modality_enums
        self.modalities = []
        self.mod_options = mod_options

    def compute_features_matrices(self, data: MultimodalData, fit: bool=True) -> List[np.ndarray]:
        matrices = []
        if fit:
            self.modalities = []
            for modality_enum in self.modality_enums:
                modality = MODALITY_CLASSES[modality_enum](use_kernel=False, **self.mod_options)
                matrices.append(modality.fit_transform(data))
                self.modalities.append(modality)
        else:
            for modality in self.modalities:
                matrices.append(modality.transform(data))
        return matrices

    def fit(self, data: MultimodalData) -> 'FeaturesExtractor':
        return self

    def fit_transform(self, data: MultimodalData) -> List[np.ndarray]:
        self.fit(data)
        return self.compute_features_matrices(data, fit=True)

    def transform(self, data: MultimodalData) -> List[np.ndarray]:
        return self.compute_features_matrices(data, fit=False)


def moments_to_beta_params(m1, m2):
    mean = m1
    var = m2 - np.square(mean)
    nu = mean * (1 - mean) / var - 1
    a = mean * nu
    b = (1 - mean) * nu
    return a, b


def dd_distance_block(X, embeddings):
    n = len(X)
    X = RobustScaler().fit_transform(X)
    embeddings = np.concatenate((embeddings, -embeddings), axis=0)
    X = np.abs(np.concatenate((np.clip(X, 0, None), np.clip(X, None, 0)), axis=1))
    X = X ** 2
    X = X / np.sum(X, axis=1)[:, np.newaxis]
    D = np.zeros((n, n))
    m = len(embeddings)
    for i in tqdm.tqdm(range(n)):

        for j in range(i):
            distance = ot.gaussian.empirical_bures_wasserstein_distance(embeddings, embeddings, ws=X[i, :, np.newaxis], wt=X[j, :, np.newaxis])
            #gamma = ot.emd(X[i, :], X[j, :], C)
            #distance = np.sum(gamma * C)
            D[i, j] = distance
            D[j, i] = D[i, j]

    D = D / np.mean(D)
    return D


def dd_distance_block(X, embeddings):
    n = len(X)
    X = RobustScaler().fit_transform(X)
    X = np.abs(X)
    X = X ** 2
    X = X / np.sum(X, axis=1)[:, np.newaxis]
    D = np.zeros((n, n))
    m = len(embeddings)
    for i in tqdm.tqdm(range(n)):

        for j in range(i):
            distance = ot.gaussian.empirical_bures_wasserstein_distance(embeddings, embeddings, ws=X[i, :, np.newaxis], wt=X[j, :, np.newaxis])
            #gamma = ot.emd(X[i, :], X[j, :], C)
            #distance = np.sum(gamma * C)
            D[i, j] = distance
            D[j, i] = D[i, j]

    D = D / np.mean(D)
    return D


def dd_distance(X):
    n = len(X)
    from sklearn.decomposition import PCA
    model = PCA(n_components=10)
    embeddings = model.fit_transform(X.T)
    embeddings = RobustScaler().fit_transform(embeddings)
    D = dd_distance_block(X, embeddings)
    return D


# https://www.sciencedirect.com/science/article/abs/pii/S0167715209001473
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224338


"""
def compute_pairwise_beta_w_distances(m1_1, m2_1, m1_2, m2_2, n_bins=100):
    n = len(m1_1)
    m = len(m1_1)
    n_features = m1_1.shape[1]
    
    distances = np.zeros((n, m))
    
    for idx in tqdm.tqdm(np.array_split(np.arange(n_features), 1000), desc='Beta W'):

        a1, b1 = moments_to_beta_params(m1_1, m2_1)
        a2, b2 = moments_to_beta_params(m1_2, m2_2)
        
        a1 = a1[:, :, np.newaxis]
        b1 = b1[:, :, np.newaxis]
        a2 = a2[:, :, np.newaxis]
        b2 = b2[:, :, np.newaxis]

        xs = np.linspace(0, 1, n_bins)[np.newaxis, np.newaxis, :]
        
        f1 = scipy.stats.beta.cdf(xs, a1, b1) / n_bins
        f2 = scipy.stats.beta.cdf(xs, a2, b2) / n_bins
        
        f1 = f1.reshape(len(f1), -1)
        f2 = f2.reshape(len(f2), -1)

        print('A', f1.shape, f2.shape)

        distances = scipy.spatial.distance.cdist(f1, f2, metric='cityblock')
    
    distances /= n_features
    
    return distances
"""


def compute_pairwise_beta_w_distances(m1_1, m2_1, m1_2, m2_2, n_bins=100):
    n = len(m1_1)
    m = len(m1_1)
    n_features = m1_1.shape[1]
    
    distances = np.zeros((n, m))
    
    for idx in tqdm.tqdm(np.array_split(np.arange(n_features), 1000), desc='Beta W'):

        a1, b1 = moments_to_beta_params(m1_1, m2_1)
        a2, b2 = moments_to_beta_params(m1_2, m2_2)
        
        a1 = a1[:, :, np.newaxis]
        b1 = b1[:, :, np.newaxis]
        a2 = a2[:, :, np.newaxis]
        b2 = b2[:, :, np.newaxis]

        xs = np.linspace(0, 1, n_bins)[np.newaxis, np.newaxis, :]
        
        f1 = scipy.stats.beta.cdf(xs, a1, b1) / n_bins
        f2 = scipy.stats.beta.cdf(xs, a2, b2) / n_bins
        
        f1 = f1.reshape(len(f1), -1)
        f2 = f2.reshape(len(f2), -1)

        distances = scipy.spatial.distance.cdist(f1, f2, metric='cityblock')
    
    distances /= n_features
    
    return distances
