from typing import Optional

import tqdm
import numpy as np
import scipy.special
import torch


def log_beta(alpha: np.ndarray) -> np.ndarray:
    alpha0 = float(np.sum(alpha))
    return np.sum(scipy.special.loggamma(alpha)) - scipy.special.loggamma(alpha0)


def inv_digamma(y: np.ndarray, x0: np.ndarray) -> np.ndarray:
    f = lambda x: scipy.special.digamma(x) - y
    f_prime = lambda x: scipy.special.polygamma(1, x)
    return scipy.optimize.newton(f, x0, fprime=f_prime)


def infer_concentration_params(X: np.ndarray) -> np.ndarray:
    X = X / np.sum(X, axis=1)[:, np.newaxis]
    mean = np.mean(X, axis=0)
    var = np.var(X, ddof=1)
    alpha_tilde = mean
    alpha_0 = (alpha_tilde * (1 - alpha_tilde)) / var - 1
    return alpha_tilde * alpha_0


def dirichlet_ll(p: np.ndarray, alpha: np.ndarray) -> float:
    return np.sum((alpha[np.newaxis, :] - 1) * np.log(p), axis=1) - log_beta(alpha)


def multinomial_ll(counts: np.ndarray, p: np.ndarray) -> float:
    return np.sum(counts * np.log(p), axis=1)


def dirichlet_multinomial_ll(counts: np.ndarray, p: np.ndarray, alpha: np.ndarray) -> float:
    return np.mean(dirichlet_ll(p, alpha) + multinomial_ll(counts, p))


def update_frequencies(counts: np.ndarray, alpha: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    a = alpha[np.newaxis, :] + counts
    a0 = np.sum(a, axis=1)[:, np.newaxis]
    p = (a - 1) / (a0 - a.shape[1])
    p = np.clip(p, eps, 1 - eps)
    p /= np.sum(p, axis=1)[:, np.newaxis]
    return p


def update_concentration_parameters(alpha: np.ndarray, p: np.ndarray) -> np.ndarray:
    # Proved in "The Book of Statistical Proofs", Soch, Joram, et al. (2024).
    log_mean = np.mean(np.log(p), axis=0)
    assert not np.any(np.isnan(log_mean))
    psi = scipy.special.digamma(np.sum(alpha))
    assert not np.any(np.isnan(psi))
    alpha = inv_digamma(psi + log_mean, alpha)
    alpha = np.clip(alpha, 1e-5, None)
    return alpha


def infer_histogram_frequencies(counts: np.ndarray, tol: float = 1e-9, eps: float = 1e-10, max_n_iter: int = 500) -> np.ndarray:

    p = counts / np.sum(counts, axis=1)[:, np.newaxis]
    p = np.clip(p, eps, 1 - eps)
    p /= np.sum(p, axis=1)[:, np.newaxis]

    alpha = infer_concentration_params(p)

    losses = []
    losses.append(-dirichlet_multinomial_ll(counts, p, alpha))

    pbar = tqdm.tqdm(range(max_n_iter), desc='Infer concentration params')
    for _ in pbar:

        p = update_frequencies(counts, alpha, eps=eps)
        alpha = update_concentration_parameters(alpha, p)
        loss = -dirichlet_multinomial_ll(counts, p, alpha)
        pbar.set_description(f'Infer concentration params (NLL={loss})')

        # Stopping criterion
        losses.append(loss)
        if abs((losses[-2] - losses[-1]) / losses[-1]) < tol:
            break

    return p, alpha



# ------------------------------------------------------------------------------------



def sigma(p_logit: np.ndarray) -> np.ndarray:
    return scipy.special.softmax(p_logit, axis=1)


def sigma_inv(p: np.ndarray) -> np.ndarray:
    return scipy.special.logit(p)


def solve_logistic_regressions_wrt_w(
        X: np.ndarray,
        p: np.ndarray,
        y: np.ndarray,
        w0: np.ndarray,
        gamma: float = 1e-1
) -> np.ndarray:

    n = np.sum(X, axis=1)[:, np.newaxis]
    k = X

    p_logit = sigma_inv(p)

    def f(w):
        pi = sigma(p_logit + y * w[np.newaxis, :])
        #return np.mean(y * (k - n * pi), axis=0) - 2 * gamma * w
        return np.sum(np.mean(y * (k - n * pi), axis=0) - 2 * gamma * w)

    def fprime(w):
        pi = sigma(p_logit + y * w[np.newaxis, :])
        return np.mean(-np.square(y) * n * pi * (1 - pi), axis=0) - 2 * gamma

    #w = scipy.optimize.newton(f, w0, fprime=fprime, maxiter=500, full_output=False, disp=True)
    w = scipy.optimize.minimize(f, w0, jac=fprime).x

    #if f(w) < f(w0):
    ##    w = w0
    #    print('WARNING')

    return w


def solve_bayesian_logistic_regressions_wrt_p(
        X: np.ndarray,
        p0: np.ndarray,
        y: np.ndarray,
        w: np.ndarray,
        alpha: np.ndarray,
        max_n_iter: int = 10,
        eps: float = 1e-10
) -> np.ndarray:

    k = X
    n = np.sum(X, axis=1)[:, np.newaxis]
    a = alpha
    b = np.sum(alpha) - alpha
    a = a[np.newaxis, :]
    b = b[np.newaxis, :]

    p0 = np.clip(p0, eps, 1 - eps)
    p0 = p0 / np.sum(p0, axis=1)[:, np.newaxis]

    p_logit = sigma_inv(p0)
    pi = sigma(p_logit + y * w[np.newaxis, :])
    pi = np.clip(pi, eps, 1 - eps)
    pi = pi / np.sum(pi, axis=1)[:, np.newaxis]

    losses = []
    for _ in range(max_n_iter):

        # Compute gradient of the logarithm of the likelihood function + Beta prior + Jacobian
        p_logit_grad = k * (1 - pi) - (n - k) * pi
        p = sigma(p_logit)
        p = np.clip(p, eps, 1 - eps)
        p = p / np.sum(p, axis=1)[:, np.newaxis]
        p_logit_grad += (a - 1) * (1 - p) - (b - 1) * p

        # Update parameters
        lr = 1e-6
        p_logit += lr * p_logit_grad

        # Deal with numerical instability
        p = sigma(p_logit)
        p = np.clip(p, eps, 1 - eps)
        p_logit = sigma_inv(p)

        # Logistic regression predictions
        pi = sigma(p_logit + y * w[np.newaxis, :])
        pi = np.clip(pi, eps, 1 - eps)
        pi = np.clip(pi, eps, 1 - eps)

        # Compute loss
        loss = -np.mean(dirichlet_ll(p, alpha) + multinomial_ll(X, pi))
        losses.append(loss)

    p = sigma(p_logit)

    return p


def log_beta_pt(alpha: torch.Tensor) -> torch.Tensor:
    alpha0 = torch.sum(alpha)
    return torch.sum(torch.lgamma(alpha)) - torch.lgamma(alpha0)


def dirichlet_ll_pt(p: torch.Tensor, alpha: torch.Tensor) -> torch.Tensor:
    return torch.sum((alpha.unsqueeze(0) - 1) * torch.log(p), dim=1) - log_beta_pt(alpha)


def multinomial_ll_pt(counts: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
    return torch.sum(counts * torch.log(p), dim=1)


def infer_histogram_frequencies_with_confounders(counts: np.ndarray, Y: np.ndarray, gamma: float = 0.1, tol: float = 1e-9, eps: float = 1e-10, max_n_iter: int = 500, lr: float = 0.001) -> np.ndarray:

    p, alpha = infer_histogram_frequencies(counts, tol=tol, eps=eps, max_n_iter=100)
    p = np.clip(p, eps, 1 - eps)
    p /= np.sum(p, axis=1)[:, np.newaxis]
    pi = np.copy(p)
    w = np.zeros(counts.shape[1])

    Y = torch.tensor(Y)
    counts = torch.tensor(counts)
    alpha_log = torch.nn.Parameter(torch.tensor(np.log(alpha)))
    p_log = torch.nn.Parameter(torch.tensor(np.log(p)))
    w = torch.nn.Parameter(torch.tensor(w))

    optimizer = torch.optim.Adam([alpha_log, p_log, w], lr=lr)

    p = torch.softmax(p_log, dim=1)
    pi = p
    alpha = torch.exp(alpha_log)

    losses = []
    losses.append((-torch.mean(dirichlet_ll_pt(p, alpha) + multinomial_ll_pt(counts, pi))).item())

    pbar = tqdm.tqdm(range(max_n_iter), desc='Infer concentration params')
    for _ in pbar:

        optimizer.zero_grad()

        p = torch.softmax(p_log, dim=1)
        pi = torch.softmax(p_log + w.unsqueeze(0) * Y, dim=1)
        alpha = torch.exp(alpha_log)

        loss = -torch.mean(dirichlet_ll_pt(p, alpha) + multinomial_ll_pt(counts, pi)) + gamma * torch.sum(torch.square(w))
        loss.backward()

        optimizer.step()

        pbar.set_description(f'Infer concentration params (NLL={loss.item()})')

        # Stopping criterion
        losses.append(loss.item())
        #if abs((losses[-2] - losses[-1]) / losses[-1]) < tol:
        #    break

    p = torch.softmax(p_log, dim=1).cpu().data.numpy()
    alpha = torch.exp(alpha_log).cpu().data.numpy()
    w = w.cpu().data.numpy()

    return p, alpha, w


class DirichletMultinomialTransformer:

    def __init__(self, tol: float = 1e-6, eps: float = 1e-10, max_n_iter: int = 500, gamma: float = 0.1, lr: float = 0.05):
        self.tol: float = tol
        self.eps: float = eps
        self.lr: float = lr
        self.max_n_iter: int = max_n_iter
        self.gamma: float = gamma
        self.alpha: Optional[np.ndarray] = None
        self.w: Optional[np.ndarray] = None

    def fit(self, counts: np.ndarray, y: Optional[np.ndarray] = None) -> 'DirichletMultinomialTransformer':
        if y is not None:
            p, alpha, w = infer_histogram_frequencies_with_confounders(counts, y, gamma=self.gamma, tol=self.tol, eps=self.eps, max_n_iter=self.max_n_iter, lr=self.lr)
            self.w = w
        else:
            p, alpha = infer_histogram_frequencies(counts, tol=self.tol, eps=self.eps, max_n_iter=self.max_n_iter)
        self.alpha = alpha
        return self

    def fit_transform(self, counts: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        if y is not None:
            p, alpha, w = infer_histogram_frequencies_with_confounders(counts, y, gamma=self.gamma, tol=self.tol, eps=self.eps, max_n_iter=self.max_n_iter, lr=self.lr)
        else:
            p, alpha = infer_histogram_frequencies(counts, tol=self.tol, eps=self.eps, max_n_iter=self.max_n_iter)
        self.alpha = alpha
        return p

    def transform(self, counts: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        assert False
        if y is not None:
            return update_frequencies_with_confounder(counts, self.alpha, y, self.w, eps=self.eps)
        else:
            return update_frequencies(counts, self.alpha, eps=self.eps)
