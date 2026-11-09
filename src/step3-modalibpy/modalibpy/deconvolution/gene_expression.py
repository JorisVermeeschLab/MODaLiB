from typing import Tuple, List, Optional

import numpy as np
import scipy.optimize
import torch
import torch.nn.functional
import tqdm
from scipy.special import logit
from sklearn.preprocessing import RobustScaler, quantile_transform

from .optimizer import Optimizer


def pearson_corrcoef(x: torch.Tensor, y: torch.Tensor, eps: float = 1e-6) -> torch.Tensor:

    #from fast_soft_sort.pytorch_ops import soft_rank, soft_sort
    #x = soft_rank(x)
    #y = soft_rank(y)

    vx = x - torch.mean(x, dim=1).unsqueeze(1)
    vy = y - torch.mean(y, dim=1).unsqueeze(1)
    cov = torch.mean(vx * vy, dim=1)
    var_x = torch.mean(torch.square(vx), dim=1) + eps
    var_y = torch.mean(torch.square(vy), dim=1) + eps
    rho = cov / torch.sqrt(var_x * var_y)
    return rho


def crossdomain_deconvolution(X: np.ndarray, R: np.ndarray, max_n_iter: int = 100) -> Tuple[np.ndarray, np.ndarray]:

    X = quantile_transform(X.T, n_quantiles=1000, output_distribution='uniform').T
    R = quantile_transform(R.T, n_quantiles=1000, output_distribution='uniform').T

    assert not np.any(np.isnan(X))
    assert not (np.any(np.isnan(R)) or np.any(np.isinf(R)))

    n = len(X)
    m = len(R)
    p = len(X)

    # Random initial solution
    alpha = np.random.rand(n, m)
    alpha /= np.sum(alpha, axis=1)[:, np.newaxis]
    
    # Convert to tensors
    alpha_logit = torch.nn.Parameter(torch.FloatTensor(np.log(alpha)))
    X = torch.FloatTensor(X)
    R = torch.FloatTensor(R)

    optimizer = Optimizer(verbose=False)
    optimizer.add([alpha_logit], 1e-3)

    patience = 30
    n_steps_without_improvement = 0
    best_loss = np.inf
    pbar = tqdm.tqdm(range(max_n_iter))
    for iteration in pbar:

        optimizer.zero_grad()

        alpha = torch.softmax(alpha_logit, dim=1)

        X_hat = torch.mm(alpha, R)

        loss = -torch.mean(torch.abs(pearson_corrcoef(X, X_hat)))
        loss.backward()
        pbar.set_description(str(loss.item()))

        # Update parameters
        optimizer.step(loss.item())

        # scheduler.step()
        if loss.item() >= best_loss:
            n_steps_without_improvement += 1
            if n_steps_without_improvement >= patience:
                break
        else:
            best_loss = loss.item()
            n_steps_without_improvement = 0

    alpha = torch.softmax(alpha_logit, dim=1).cpu().data.numpy()

    return alpha
