from typing import Optional, Tuple

import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import BetaBinomialTransformer
from ..end_motifs import cgnngc_ratios, semfar_scores


class JaggednessModality(BaseModality):

    # Minimum fragment length for which jag index will be estimated
    LB = 70

    # Maximum fragment length for which jag index will be estimated
    UB = 200

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='rbf', **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.bb_transformer1: BetaBinomialTransformer = BetaBinomialTransformer(gamma=0.01, max_n_iter=2000, lr=1e-6)
        self.a, self.b = [], []

    @staticmethod
    def data_imputation(X: np.ndarray) -> None:
        mu = np.nanmean(X, axis=0)
        idx = np.where(~np.isnan(mu))[0]
        default_mu = mu[idx[0]]
        mask = np.isnan(X)
        mu_mask = np.isnan(mu)
        for j in range(X.shape[1]):
            for i in range(X.shape[0]):
                if mask[i, j]:
                    if j > 0:
                        X[i, j] = X[i, j - 1]
                    elif not mu_mask[j]:
                        X[i, j] = mu[j]
                    else:
                        X[i, j] = default_mu

    @staticmethod
    def extract_meth_per_length(data: MultimodalData) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        eps = 0.01
        M = data["Overall-methylated-counts"]
        U = data["Overall-unmethylated-counts"]
        T = M + U
        R = M / T
        s = 30
        X1, X2, N1, N2 = [], [], [], []
        for fragment_length in range(0, 501):
            X1.append(np.mean(R[:, fragment_length, :min(fragment_length, s)], axis=1))
            X2.append(np.mean(R[:, fragment_length, max(0, fragment_length-s):fragment_length], axis=1))
            N1.append(np.sum(T[:, fragment_length, :min(fragment_length, s)], axis=1))
            N2.append(np.sum(T[:, fragment_length, max(0, fragment_length-s):fragment_length], axis=1))
        X1 = np.asarray(X1).T
        X2 = np.asarray(X2).T
        N1 = np.asarray(N1).T
        N2 = np.asarray(N2).T

        # Data imputation
        JaggednessModality.data_imputation(X1)
        JaggednessModality.data_imputation(X2)

        return X1, X2, N1, N2

    @staticmethod
    def compute_jaggedness(data: MultimodalData, method: str = 'corrected', probabilistic: bool = False) -> np.ndarray:

        y_cgn_ncg = cgnngc_ratios(data['End-motif probability (5 prime)'])

        r1_meth = 0.5 * (np.mean(data['OT-R1'][:, :30], axis=1) + np.mean(data['OB-R1'][:, -30:], axis=1))
        r2_meth = 0.5 * (np.mean(data['OT-R2'][:, :-30], axis=1) + np.mean(data['OB-R2'][:, :30], axis=1))

        #assert not np.any(r1_meth == 0)
        r1_meth[r1_meth == 0] = np.mean(r1_meth != 0)
        r2_meth[r2_meth == 0] = np.mean(r2_meth != 0)

        jag = 100 * (r1_meth - r2_meth) / r1_meth

        # Correct for SEMFaR score
        if method == 'corrected':
            X = data['End-motif probability (5 prime)']
            X = X / np.sum(X, axis=1)[:, np.newaxis]
            c = semfar_scores(X)
            jag_hat = -30517.79899074 * c + 165.4546993
            jag = jag / jag_hat

        return jag

    @property
    def key1(self) -> str:
        return 'Methylation R1 (preprocessed)'

    @property
    def key2(self) -> str:
        return 'Methylation R2 (preprocessed)'

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if (self.key1 in data) and (self.key2 in data) and self.cache:
            return None

        X1, X2, N1, N2 = JaggednessModality.extract_meth_per_length(data)

        lb, ub = JaggednessModality.LB, JaggednessModality.UB

        y_cgn_ncg = cgnngc_ratios(data['End-motif probability (5 prime)'])
        Y = np.repeat(y_cgn_ncg[:, np.newaxis], X1.shape[1], axis=1)
        Y = Y[:, lb:ub]

        if self.probabilistic:
            if fit:
                X1[:, lb:ub] = self.bb_transformer1.fit_transform(N1[:, lb:ub], X1[:, lb:ub] * N1[:, lb:ub], y=Y)
                X2[:, lb:ub] = self.bb_transformer1.transform(N2[:, lb:ub], X2[:, lb:ub] * N2[:, lb:ub], y=Y)
            else:
                X1[:, lb:ub] = self.bb_transformer1.transform(N1[:, lb:ub], X1[:, lb:ub] * N1[:, lb:ub], y=Y)
                X2[:, lb:ub] = self.bb_transformer1.transform(N2[:, lb:ub], X2[:, lb:ub] * N2[:, lb:ub], y=Y)

        self.save(data, self.key1, X1)
        self.save(data, self.key2, X2)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X1 = self.load(data, self.key1)
        X2 = self.load(data, self.key2)

        # Compute jag index for each fragment length separately
        lb, ub = JaggednessModality.LB, JaggednessModality.UB
        r1_meth = X1[:, lb:ub]
        r2_meth = X2[:, lb:ub]
        jag = 100.0 * (r1_meth - r2_meth) / r1_meth
        jag_corrected = np.copy(jag)

        # Correct for SEMFaR score
        X = data['End-motif probability (5 prime)']
        X = X / np.sum(X, axis=1)[:, np.newaxis]
        c = semfar_scores(X)
        if fit:
            self.a, self.b = [], []
        for k in range(r1_meth.shape[1]):
            if fit:
                if (y is not None):
                    coef = np.polyfit(c[y == 0], jag[y == 0, k], 1)
                    self.a.append(coef[0])
                    self.b.append(coef[1])
                else:
                    self.a.append(-30517.79899074)
                    self.b.append(165.4546993)
            jag_hat = self.a[k] * c + self.b[k]
            jag_corrected[:, k] /= jag_hat

        # Use both jag index and S-jag index
        jag = np.concatenate((jag, jag_corrected), axis=1)

        return jag
