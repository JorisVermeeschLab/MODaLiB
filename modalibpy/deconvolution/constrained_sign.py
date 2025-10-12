from typing import Tuple, Optional

import numpy as np
import tqdm
import numpy as np
import scipy.optimize
import scipy.special
from scipy.optimize import nnls
from scipy.special import logit
import torch
import torch.nn.functional
from sklearn.decomposition import NMF

from modalibpy.deconvolution.optimizer import Optimizer


def pearson_corrcoef(x: torch.Tensor, y: torch.Tensor, eps: float = 1e-6) -> torch.Tensor:
    vx = x - torch.mean(x, dim=1).unsqueeze(1)
    vy = y - torch.mean(y, dim=1).unsqueeze(1)
    cov = torch.mean(vx * vy, dim=1)
    var_x = torch.mean(torch.square(vx), dim=1) + eps
    var_y = torch.mean(torch.square(vy), dim=1) + eps
    rho = cov / torch.sqrt(var_x * var_y)
    return rho


def constrained_sign_deconvolution(
        X: np.ndarray,
        R_sign: np.ndarray,
        max_n_iter: int = 5000,
        R: Optional[np.ndarray] = None,
        alpha_mean: Optional[np.ndarray] = None,
        gamma: float = 0.1
) -> Tuple[np.ndarray, np.ndarray]:

    # Initialize factor matrices
    fit_R = (R is None)
    if R is None:
        R = np.abs(np.random.normal(0, 1, size=R_sign.shape)) * R_sign
        alpha = np.random.rand(X.shape[0], R_sign.shape[0])
    else:
        alpha = []
        for i in range(len(X)):
            alpha.append(scipy.optimize.nnls(R.T, X[i, :])[0])
        alpha = np.asarray(alpha)
    alpha = alpha + 1e-15
    alpha_logit = torch.nn.Parameter(torch.FloatTensor(np.log(alpha / np.sum(alpha, axis=1)[:, np.newaxis])))
    R = torch.nn.Parameter(torch.FloatTensor(R))
    R_sign = torch.FloatTensor(R_sign.astype(float))

    X = torch.FloatTensor(X)
    if alpha_mean is not None:
        alpha_mean = torch.FloatTensor(alpha_mean).unsqueeze(0)

    optimizer = Optimizer(verbose=False)
    optimizer.add([alpha_logit], 1e-2)
    if fit_R:
        optimizer.add([R], 1e-2)

    patience = 1000
    n_steps_without_improvement = 0
    best_loss = np.inf
    pbar = tqdm.tqdm(range(max_n_iter))
    for iteration in pbar:

        optimizer.zero_grad()

        alpha = torch.softmax(alpha_logit, dim=1)

        R_signed = R_sign * torch.abs(R)
        X_hat = torch.mm(alpha, R_signed)

        loss = -torch.mean(pearson_corrcoef(X, X_hat))
        if alpha_mean is not None:
            loss = loss + gamma * torch.mean(torch.square(alpha - alpha_mean))
        loss = loss + gamma * torch.mean(torch.square(torch.clamp(0.03 - alpha, 0, None)))

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
    R = torch.abs(R).cpu().data.numpy()

    return alpha, R


def unconstrained_sign_deconvolution(
        X: np.ndarray,
        y: np.ndarray,
        max_n_iter: int = 1000,
) -> Tuple[np.ndarray, np.ndarray]:

    # Initialize factor matrices
    alpha = torch.nn.Parameter(torch.FloatTensor(np.random.normal(0, 1, size=(X.shape[1], y.shape[1]))))
    print(alpha.size())

    # Initialize bias term
    b = np.full(X.shape[1], 1e-5)
    bias = torch.nn.Parameter(torch.FloatTensor(b))

    X = torch.FloatTensor(X)
    y = torch.FloatTensor(y)

    optimizer = Optimizer(verbose=False)
    optimizer.add([alpha], 1e-2)
    optimizer.add([bias], 1e-2)

    patience = 100
    n_steps_without_improvement = 0
    best_loss = np.inf
    pbar = tqdm.tqdm(range(max_n_iter))
    for iteration in pbar:

        optimizer.zero_grad()

        y_hat = torch.mm(X + bias.unsqueeze(0), alpha)
        loss = torch.mean(torch.square(y - y_hat))
        loss.backward()
        pbar.set_description(str(loss))

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

    alpha = torch.abs(alpha).cpu().data.numpy()
    b = bias.cpu().data.numpy()

    return alpha, b
