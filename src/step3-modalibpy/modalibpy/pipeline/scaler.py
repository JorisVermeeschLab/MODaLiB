import numpy as np
import scipy.spatial
from sklearn.preprocessing import StandardScaler, RobustScaler

from modalibpy.pipeline.base import BaseStep


class Scaler(BaseStep):

    def __init__(self):
        BaseStep.__init__(self)
        self.scaler = RobustScaler()

    def fit(self, X):
        self.scaler.fit(X)
        return self

    def transform(self, dist, X):
        X = self.scaler.transform(X)

        dist_transformed = dist - self.scaler.center_[np.newaxis, :]
        dist_transformed = dist_transformed * (1. / self.scaler.scale_[np.newaxis, :])

        if not self.compare(X, dist, dist_transformed):
            dist_transformed = dist

        return dist_transformed, X

    def fit_transform(self, dist, X):
        self.fit(X)
        return self.transform(dist, X)
