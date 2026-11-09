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

        if ('Fragment size index (preprocessed)' in data) and self.cache:
            return None

        if 'Fragment size index' in data:
            X = data['Fragment size index']
        else:
            X = data['Proportion of long fragments']
        n = data['Coverage']

        self.save(data, 'Fragment size index (preprocessed)', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = self.load(data, 'Fragment size index (preprocessed)')
        n = np.clip(data['Coverage'], 0.01, None)

        return X
