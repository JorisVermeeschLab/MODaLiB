import numpy as np
import torch
import pyro
from torch.distributions import constraints, Distribution
from pyro.distributions import BetaBinomial

from modalibpy.inference.svi import svi_optimize


def __model_mle(n: torch.Tensor, k: torch.Tensor) -> None:

    # Ensure k is in the support
    k = torch.round(k)

    # Initial parameter estimation
    p = (k / n).cpu().data.numpy()
    mean = np.mean(p, axis=0)
    var = np.var(p, axis=0, ddof=1)
    shared = mean * (1. - mean) / var - 1
    a_hat = mean * shared
    b_hat = (1. - mean) * shared
    assert not (np.any(np.isnan(a_hat)) or np.any(np.isinf(a_hat)))
    assert not (np.any(np.isnan(b_hat)) or np.any(np.isinf(b_hat)))

    a = pyro.param('alpha', torch.tensor(a_hat), constraint=constraints.positive)
    b = pyro.param('beta', torch.tensor(b_hat), constraint=constraints.positive)
    with pyro.plate('data', n.size(1)):
        dist = BetaBinomial(a, b, total_count=n)
        pyro.sample('obs', dist, obs=k)


def infer_beta_binomial(n: np.ndarray, k: np.ndarray) -> Distribution:
    svi_optimize([n, k], __model_mle, n_steps=30, lr=0.00002)  # TODO
    a = pyro.param('alpha')
    b = pyro.param('beta')
    dist = BetaBinomial(a.data, b.data, total_count=torch.FloatTensor(n))
    return dist
