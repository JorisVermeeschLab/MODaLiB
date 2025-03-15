from typing import Optional

import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import BetaBinomialTransformer
from ..end_motifs import cgnngc_ratios


class JaggednessModality(BaseModality):

    # Minimum fragment length for which jag index will be estimated
    LB = 330

    # Maximum fragment length for which jag index will be estimated
    UB = 380

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='rbf', **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.bb_transformer1: BetaBinomialTransformer = BetaBinomialTransformer(gamma=0.01, max_n_iter=2000, lr=1e-6)

    @staticmethod
    def compute_jaggedness(data: MultimodalData, method: str = 'corrected', probabilistic: bool = False) -> np.ndarray:

        y_cgn_ncg = cgnngc_ratios(data['End-motif probability (5 prime)'])

        if method == 'jag-index':
            r1_meth = 0.5 * (np.mean(data['OT-R1'][:, :30], axis=1) + np.mean(data['OB-R1'][:, -30:], axis=1))
            r2_meth = 0.5 * (np.mean(data['OT-R2'][:, :-30], axis=1) + np.mean(data['OB-R2'][:, :30], axis=1))
        elif method == 'corrected':
            lb, ub = JaggednessModality.LB, JaggednessModality.UB
            X1 = data['R1 methylation'][:, lb:ub]
            X2 = data['R2 methylation'][:, lb:ub]
            if probabilistic:
                n = data['Fragment length counts'][:, lb:ub] * 10
                Y = np.repeat(y_cgn_ncg[:, np.newaxis], X.shape[1], axis=1)
                bb_transformer = BetaBinomialTransformer(gamma=0.01)
                X1 = bb_transformer.fit_transform(n, n * X1, y=Y)
                X2 = bb_transformer.transform(n, n * X2, y=Y)
            r1_meth = np.mean(X1, axis=1)
            r2_meth = np.mean(X2, axis=1)
        else:
            raise NotImplementedError(f'Unknown jaggedness estimation method "{method}"')
        #assert not np.any(r1_meth == 0)
        r1_meth[r1_meth == 0] = np.mean(r1_meth != 0)
        r2_meth[r2_meth == 0] = np.mean(r2_meth != 0)

        jag = (r1_meth - r2_meth) / r1_meth

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

        X1 = np.copy(data['R1 methylation'])
        X2 = np.copy(data['R2 methylation'])
        n = np.clip(data['Fragment length counts'], 0.1, None) * 5

        lb, ub = JaggednessModality.LB, JaggednessModality.UB

        y_cgn_ncg = cgnngc_ratios(data['End-motif probability (5 prime)'])
        Y = np.repeat(y_cgn_ncg[:, np.newaxis], X1.shape[1], axis=1)
        Y = Y[:, lb:ub]

        if self.probabilistic:
            if fit:
                X1[:, lb:ub] = self.bb_transformer1.fit_transform(n[:, lb:ub], X1[:, lb:ub] * n[:, lb:ub], y=Y)
                X2[:, lb:ub] = self.bb_transformer1.transform(n[:, lb:ub], X2[:, lb:ub] * n[:, lb:ub], y=Y)
            else:
                X1[:, lb:ub] = self.bb_transformer1.transform(n[:, lb:ub], X1[:, lb:ub] * n[:, lb:ub], y=Y)
                X2[:, lb:ub] = self.bb_transformer1.transform(n[:, lb:ub], X2[:, lb:ub] * n[:, lb:ub], y=Y)

        self.save(data, self.key1, X1)
        self.save(data, self.key2, X2)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X1 = self.load(data, self.key1)
        X2 = self.load(data, self.key2)

        # Compute jag index for each fragment length separately
        lb, ub = JaggednessModality.LB, JaggednessModality.UB
        r1_meth = X1[:, lb:ub]
        r2_meth = X2[:, lb:ub]
        X = (r1_meth - r2_meth) / r1_meth

        X = np.nan_to_num(X, nan=0, posinf=0, neginf=0)

        return X
