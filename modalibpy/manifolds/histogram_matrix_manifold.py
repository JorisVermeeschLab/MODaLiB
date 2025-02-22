import torch

from modalibpy.manifolds.base import MatrixManifold


class HistogramMatrixManifold(MatrixManifold):

    def __init__(self, *args, **kwargs):
        MatrixManifold.__init__(self, *args, **kwargs)

    def initialize(self) -> torch.Tensor:
        R = torch.rand(self.n, self.m)
        R = R / torch.sum(R, dim=1).unsqueeze(1)
        return torch.log(R)

    def retraction(self, R: torch.Tensor) -> torch.Tensor:
        return torch.softmax(R, dim=1)

    @torch.no_grad()
    def _contains(self, R: torch.Tensor) -> bool:
        if not torch.all(R >= -self.eps):
            return False
        if not torch.all(torch.abs(torch.sum(R, dim=1) - 1) <= eps):
            return False
        return True
