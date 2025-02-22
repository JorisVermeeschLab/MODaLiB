import numpy as np

from modalibpy.distributions.base import BaseDistribution


class GaussianDistribution(BaseDistribution):

    def __init__(self, mu: np.ndarray, sigma2: np.ndarray):
        self.mu: np.ndarray
        self.sigma2: np.ndarray

    def mean(self) -> np.ndarray:
        return self.mu

    def variance(self) -> np.ndarray:
        return self.sigma2
