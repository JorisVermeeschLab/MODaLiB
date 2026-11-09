import torch

from modalibpy.manifolds.base import MatrixManifold


class RatioMatrixManifold(MatrixManifold):

    def __init__(self, *args, **kwargs):
        MatrixManifold.__init__(self, *args, **kwargs)

    def initialize(self) -> torch.Tensor:
        R = torch.rand(self.n, self.m)
        return torch.logit(R)

    def retraction(self, R: torch.Tensor) -> torch.Tensor:
        return torch.sigmoid(R)

    @torch.no_grad()
    def _contains(self, R: torch.Tensor) -> bool:
        if not torch.all(R >= -self.eps):
            return False
        if not torch.all(R <= 1 + self.eps):
            return False
        return True
