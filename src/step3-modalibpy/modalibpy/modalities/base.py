from typing import Dict, Optional
from abc import ABCMeta, abstractmethod

import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import KernelCenterer

from ..kernels.rbf import RBFKernel
from ..kernels.hmm import HMMKernel
from ..kernels.chi_squared import ChiSquaredKernel
from ..kernels.wasserstein import WassersteinKernel, GaussianWassersteinKernel
from ..kernels.anova import *
from ..kernels.braycurtis import *


MultimodalData = Dict[str, np.ndarray]


class BaseModality(object):

    def __init__(self, use_kernel: bool = True, kernel: str = 'rbf', probabilistic: bool = False, cache: bool = False, **kwargs) -> None:
        self.kernel_type: str = kernel
        if self.kernel_type == 'rbf':
            self.kernel = RBFKernel()
        elif self.kernel_type == 'chi2':
            self.kernel = ChiSquaredKernel()
        elif self.kernel_type == 'hmm':
            self.kernel = HMMKernel()
        elif self.kernel_type == 'anova':
            self.kernel = ANOVAKernel()
        elif self.kernel_type == 'wasserstein':
            self.kernel = WassersteinKernel()
        elif self.kernel_type == 'gaussian-wasserstein':
            self.kernel = GaussianWassersteinKernel()
        elif self.kernel_type == 'braycurtis':
            self.kernel = BrayCurtisKernel()
        else:
            raise NotImplementedError(f'Unknown kernel "{self.kernel_type}"')
        self.use_kernel: bool = use_kernel

        self.probabilistic: bool = probabilistic
        self.cache: bool = cache

    def fit(self, data: MultimodalData, y: Optional[np.ndarray] = None) -> 'BaseModality':
        X = self.process(data, y=y, fit=True)
        if self.use_kernel:
            self.kernel.fit(X)
        return self

    def transform(self, data: MultimodalData) -> np.ndarray:
        X = self.process(data, fit=False)
        if self.use_kernel:
            X = self.kernel.transform(X)
        return X

    def fit_transform(self, data: MultimodalData, y: Optional[np.ndarray] = None) -> np.ndarray:
        X = self.process(data, y=y, fit=True)
        if self.use_kernel:
            X = self.kernel.fit_transform(X)
        return X

    @abstractmethod
    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:
        pass

    @abstractmethod
    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        pass

    def load(self, data: MultimodalData, key: str) -> Optional[np.ndarray]:
        if key not in data:
            self.preprocess(data, fit=True)
        return data[key]

    def save(self, data: MultimodalData, key: str, X: np.ndarray) -> None:
        data[key] = X
