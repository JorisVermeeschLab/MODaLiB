from typing import Optional

import numpy as np
from sklearn.preprocessing import RobustScaler, StandardScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import BetaBinomialTransformer
from ..decorrelate import gc_correction
from ..end_motifs import *
from ..pipeline.fs import FeatureSelector


class LongFragmentRatioModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.bb_transformer: BetaBinomialTransformer = BetaBinomialTransformer(gamma=0.01)
        self.mask: Optional[np.ndarray] = None

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('Proportion of long fragments (preprocessed)' in data) and self.cache:
            return None

        X = data['Proportion of long fragments']
        n = data['Coverage']

        # Remove low-coverage bins
        #if fit:
        #    self.mask = (np.mean(n, axis=0) > 200)
        #X = X[:, self.mask]
        #n = n[:, self.mask]

        #if self.probabilistic:
        #    y = semfar_scores(data['End-motif probability (5 prime)'])
        #    Y = np.repeat(y[:, np.newaxis], X.shape[1], axis=1)
        #    Y = StandardScaler().fit_transform(Y)
        #    if fit:
        #        X = self.bb_transformer.fit_transform(n, X * n, y=Y)
        #    else:
        #        X = self.bb_transformer.transform(n, X * n, y=Y)

        self.save(data, 'Proportion of long fragments (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = self.load(data, 'Proportion of long fragments (preprocessed)')
        n = np.clip(data['Coverage'], 0.01, None)

        return X
