import numpy as np

from modalibpy.distributions.base import BaseDistribution


class BetaBinomialModel(BaseDistribution):

    def __init__(self, n: np.ndarray, p_mean: np.ndarray, p_variance: np.ndarray):
        alpha, beta = BetaBinomialModel.beta_parameters_from_mean_and_variance(p_mean, p_variance)
        self.n: np.ndarray = n
        self.alpha: np.ndarray = alpha
        self.beta: np.ndarray = beta

    def mean(self) -> np.ndarray:
        return self.n * self.alpha / (self.alpha + self.beta)

    def variance(self) -> np.ndarray:
        numerator = self.n * self.alpha * self.beta # * (self.alpha + self.beta + self.n)
        denominator = np.square(self.alpha + self.beta) * (self.alpha + self.beta + 1)
        return numerator / denominator

    @staticmethod
    def beta_parameters_from_mean_and_variance(mean: np.ndarray, var: np.ndarray, eps: float = 1e-20):
        shared = (mean * (1.0 - mean)) / var - 1
        alpha = np.clip(mean * shared, eps, None)
        beta = np.clip((1.0 - mean) * shared, eps, None)
        return alpha, beta
