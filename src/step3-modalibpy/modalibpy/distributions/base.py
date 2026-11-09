from abc import ABCMeta, abstractmethod

import numpy as np


class BaseDistribution(metaclass=ABCMeta):

    @abstractmethod
    def mean(self) -> np.ndarray:
        pass

    @abstractmethod
    def variance(self) -> np.ndarray:
        pass

    def std(self) -> np.ndarray:
        return np.sqrt(self.variance())

    def sqeuclidean_cdist(self) -> np.ndarray:
        mu = self.mean()
        sigma2 = self.variance()
        e2 = np.sum(sigma2 + np.square(mu), axis=1)[:, np.newaxis]
        return e2 + e2.T - 2 * np.dot(mu, mu.T)

    def __add__(self, intercept) -> 'BaseDistribution':
        return AffineTransformation(self, intercept=intercept)

    def __radd__(self, intercept) -> 'BaseDistribution':
        return self.__add__(intercept)

    def __iadd__(self, intercept) -> 'BaseDistribution':
        raise NotImplementedError()

    def __sub__(self, intercept) -> 'BaseDistribution':
        return AffineTransformation(self, intercept=-intercept)

    def __rsub__(self, intercept) -> 'BaseDistribution':
        raise NotImplementedError()

    def __isub__(self, intercept) -> 'BaseDistribution':
        raise NotImplementedError()

    def __mul__(self, scale) -> 'BaseDistribution':
        return AffineTransformation(self, scale=scale)

    def __rmul__(self, scale) -> 'BaseDistribution':
        return self.__mul__(scale)

    def __imul__(self, scale) -> 'BaseDistribution':
        raise NotImplementedError()

    def __div__(self, scale) -> 'BaseDistribution':
        return AffineTransformation(self, scale=1./scale)

    def __rdiv__(self, scale) -> 'BaseDistribution':
        return self.__div__(scale)

    def __idiv__(self, scale) -> 'BaseDistribution':
        raise NotImplementedError()


class AffineTransformation(BaseDistribution):

    def __init__(self, dist: BaseDistribution, scale=1, intercept=0):
        self.dist = dist
        self.scale = scale
        self.intercept = intercept

    def mean(self) -> np.ndarray:
        return self.intercept + self.scale * self.dist.mean()

    def variance(self) -> np.ndarray:
        return np.square(self.scale) * self.dist.variance()
