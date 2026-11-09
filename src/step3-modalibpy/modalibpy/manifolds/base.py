from abc import ABCMeta, abstractmethod

import torch


class MatrixManifold(metaclass=ABCMeta):

    def __init__(self, n: int, m: int, eps: float = 1e-7):
        self.n: int = n
        self.m: int = m
        self.eps: float = eps

    @abstractmethod
    def initialize(self) -> torch.Tensor:
        pass

    @abstractmethod
    def retraction(self, R: torch.Tensor) -> torch.Tensor:
        pass

    @torch.no_grad()
    def __contains__(self, R: torch.Tensor) -> bool:
        if len(R.size()) != 2:
            return False
        if (R.size()[0] != self.n) or (R.size()[1] != self.m):
            return False
        return self._contains(R)

    @abstractmethod
    def _contains(self, R: torch.Tensor) -> bool:
        pass
