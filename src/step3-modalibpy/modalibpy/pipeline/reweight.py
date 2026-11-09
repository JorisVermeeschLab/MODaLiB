import numpy as np
import scipy.spatial

from modalibpy.pipeline.base import BaseStep


class Reweight(BaseStep):

    def __init__(self):
        BaseStep.__init__(self)

    def fit_transform(self, dist, X, weights):

        if len(weights.shape) == 1:
            weights = weights[np.newaxis, :]

        X = X * weights
        dist_transformed = dist * weights

        if not self.compare(X, dist, dist_transformed):
            dist_transformed = dist

        return dist_transformed, X
