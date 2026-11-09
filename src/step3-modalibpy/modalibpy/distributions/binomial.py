import numpy as np

from modalibpy.distributions.base import BaseDistribution


class BinomialDistribution(BaseDistribution):

    def __init__(self, n: np.ndarray, p: np.ndarray):
        self.n: np.ndarray = n
        self.p: np.ndarray = p

    def mean(self) -> np.ndarray:
        return self.n * self.p

    def variance(self) -> np.ndarray:
        return self.n * self.p * (1. - self.p)
