import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData


class GeneExpressionModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, **kwargs)
        self.scaler: RobustScaler = RobustScaler()

    def process(self, data: MultimodalData, fit: bool = True) -> np.ndarray:
        X = data['gene-expression']
        if fit:
            X = self.scaler.fit_transform(X)
        else:
            X = self.scaler.transform(X)
        return X
