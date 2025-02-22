from typing import Optional

import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData
from ..preprocessing import infer_beta_values
from ..decorrelate import gc_correction


class CHMethylationModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.feature_weights: Optional[np.ndarray] = None

    def process(self, data: MultimodalData, fit: bool = True) -> np.ndarray:

        X = self.load(data, 'CH methylation (cache)')
        if X is not None:
            return X

        X = data['CH Methylation (first moment)']

        #X = infer_beta_values(n, k)

        # X = gc_correction(X, data['GC-content'])

        if fit:
            X = self.scaler.fit_transform(X)
        else:
            X = self.scaler.transform(X)

        self.save(data, 'CH methylation (cache)', X)

        return X
