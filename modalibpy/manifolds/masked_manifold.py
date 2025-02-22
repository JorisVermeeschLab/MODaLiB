import numpy as np
import torch

from modalibpy.manifolds.base import MatrixManifold


class MaskedManifold(MatrixManifold):

    def __init__(self, mask: np.ndarray, values: np.ndarray, manifold: MatrixManifold, **kwargs):
        MatrixManifold.__init__(self, manifold.n, manifold.m, **kwargs)
        self.mask: torch.BoolTensor = torch.BoolTensor(mask.astype(bool))
        self.values: torch.Tensor = torch.FloatTensor(values.astype(float))
        self.manifold: MatrixManifold = manifold

    def initialize(self) -> torch.Tensor:
        R = self.manifold.initialize()
        return self.mask * self.values + (~self.mask) * R

    def retraction(self, R: torch.Tensor) -> torch.Tensor:
        R = self.manifold.retraction(R)
        return self.mask * self.values + (~self.mask) * R

    @torch.no_grad()
    def _contains(self, R: torch.Tensor) -> bool:
        # TODO
        return torch.all(R[self.mask] == self.values)
