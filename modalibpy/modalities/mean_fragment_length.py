import numpy as np
from sklearn.preprocessing import RobustScaler

from .base import BaseModality, MultimodalData


class MeanFragmentLengthModality(BaseModality):

    def __init__(self, *args, **kwargs) -> None:
        BaseModality.__init__(self, *args, **kwargs)
        self.scaler: RobustScaler = RobustScaler()

    def process(self, data: MultimodalData, fit: bool = True) -> np.ndarray:

        X = self.load(data, 'Mean fragment lengths (cache)')
        if X is not None:
            return X

        X = data['Fragment length (first moment)']
        if fit:
            X = self.scaler.fit_transform(X)
        else:
            X = self.scaler.transform(X)

        self.save(data, 'Mean fragment lengths (cache)', X)

        return X
