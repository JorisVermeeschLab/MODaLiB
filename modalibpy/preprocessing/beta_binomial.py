from typing import Tuple, Optional

import tqdm
import numpy as np
import scipy.special


def beta_ll(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.newaxis, :]
    b = b[np.newaxis, :]
    return (a - 1) * np.log(p) + (b - 1) * np.log(1 - p) - scipy.special.betaln(a, b)


def binomial_ll(n: np.ndarray, k: np.ndarray, p: np.ndarray) -> float:
    return k * np.log(p) + (n - k) * np.log(1 - p)


def beta_binomial_ll(n: np.ndarray, k: np.ndarray, p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    binomial_pmf = binomial_ll(n, k, p)
    beta_pdf = beta_ll(p, a, b)
    return np.mean(binomial_pmf + beta_pdf)


def beta_parameters_from_beta_values(p: np.ndarray, eps: float = 1e-1) -> np.ndarray:
    mean = np.mean(p, axis=0)
    var = np.var(p, axis=0, ddof=1)
    var = np.clip(var, 1e-4, None)
    shared = mean * (1. - mean) / var - 1
    a = mean * shared
    b = (1. - mean) * shared
    a = np.clip(a, eps, None)
    b = np.clip(b, eps, None)
    return a, b


def mean_and_variance_from_beta_parameters(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mean = a / (a + b)
    var = (a * b) / (np.square(a + b) * (a + b + a))
    return mean, var


def beta_values_from_beta_parameters(n: np.ndarray, k: np.ndarray, a: np.ndarray, b: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    alpha = k + a
    beta = n - k + b
    p = alpha / (alpha + beta)
    p = np.clip(p, eps, 1 - eps)
    return p


def infer_beta_values(n: np.ndarray, k: np.ndarray, tol: float = 1e-6, max_n_iter: int = 100) -> np.ndarray:
    
    p = k / n
    eps = 1e-7
    p = np.clip(p, eps, 1 - eps)
    a, b = beta_parameters_from_beta_values(p)

    losses = []
    losses.append(-beta_binomial_ll(n, k, p, a, b))

    pbar = tqdm.tqdm(range(max_n_iter), desc='Infer beta values')
    for _ in pbar:

        p = beta_values_from_beta_parameters(n, k, a, b)
        a, b = beta_parameters_from_beta_values(p)
        loss = -beta_binomial_ll(n, k, p, a, b)
        pbar.set_description(f'Infer beta values (NLL={loss})')

        # Stopping criterion
        losses.append(loss)
        if abs((losses[-2] - losses[-1]) / losses[-1]) < tol:
            break

    return p, a, b


def solve_logistic_regressions_wrt_w(
        n: np.ndarray,
        k: np.ndarray,
        p: np.ndarray,
        y: np.ndarray,
        w0: np.ndarray,
        gamma: float = 1e-1
) -> np.ndarray:

    p_logit = scipy.special.logit(p)

    def f(w):
        pi = scipy.special.expit(p_logit + y * w[np.newaxis, :])
        return np.mean(y * (k - n * pi), axis=0) - 2 * gamma * w

    def fprime(w):
        pi = scipy.special.expit(p_logit + y * w[np.newaxis, :])
        return np.mean(-np.square(y) * n * pi * (1 - pi), axis=0) - 2 * gamma

    w = scipy.optimize.newton(f, w0, fprime=fprime)

    return w


def solve_bayesian_logistic_regressions_wrt_p(
        n: np.ndarray,
        k: np.ndarray,
        p0: np.ndarray,
        y: np.ndarray,
        w: np.ndarray,
        a: np.ndarray,
        b: np.ndarray,
        max_n_iter: int = 10,
        eps: float = 1e-7
) -> np.ndarray:

    a = a[np.newaxis, :]
    b = b[np.newaxis, :]

    p0 = np.clip(p0, eps, 1 - eps)
    p_logit = scipy.special.logit(p0)
    pi = scipy.special.expit(p_logit + y * w[np.newaxis, :])
    pi = np.clip(pi, eps, 1 - eps)

    losses = []
    loss = -np.mean(binomial_ll(n, k, pi) + beta_ll(p0, a, b))
    losses.append(loss)

    for _ in range(max_n_iter):

        # Compute gradient of the logarithm of the likelihood function + Beta prior + Jacobian
        p_logit_grad = k * (1 - pi) - (n - k) * pi
        p = scipy.special.expit(p_logit)
        p = np.clip(p, eps, 1 - eps)
        p_logit_grad += (a - 1) * (1 - p) - (b - 1) * p

        # Update parameters
        lr = 2 * 1e-6
        p_logit += lr * p_logit_grad

        # Deal with numerical instability
        p = scipy.special.expit(p_logit)
        p = np.clip(p, eps, 1 - eps)
        p_logit = scipy.special.logit(p)

        # Logistic regression predictions
        pi = scipy.special.expit(p_logit + y * w[np.newaxis, :])
        pi = np.clip(pi, eps, 1 - eps)

        # Compute loss
        loss = -np.mean(binomial_ll(n, k, pi) + beta_ll(p, a, b))
        losses.append(loss)

    improving = (losses[-1] < losses[0])
    # assert improving

    return p


def infer_beta_values_with_confounders(n: np.ndarray, k: np.ndarray, Y: np.ndarray, tol: float = 1e-10, max_n_iter: int = 500, gamma: float = 0.1) -> np.ndarray:

    p = k / n
    eps = 1e-7
    p, a, b = infer_beta_values(n, k, tol=tol, max_n_iter=max_n_iter)
    pi = np.copy(p)
    w = np.zeros(n.shape[1])

    losses = []
    losses.append(np.mean(-beta_ll(p, a, b) - binomial_ll(n, k, pi)))

    pbar = tqdm.tqdm(range(max_n_iter), desc='Infer beta values')
    for _ in pbar:

        # Solve Bayesian logistic regressions w.r.t. p
        p_new = solve_bayesian_logistic_regressions_wrt_p(n, k, p, Y, w, a, b, eps=eps)
        p_new = np.clip(p_new, eps, 1 - eps)

        # Solve logistic regressions w.r.t. w
        p_logit = scipy.special.logit(p)
        w_new = solve_logistic_regressions_wrt_w(n, k, p, Y, w, gamma=gamma)
        pi_new = np.clip(scipy.special.expit(p_logit + Y * w[np.newaxis, :]), eps, 1 - eps)

        # Update latent variables
        p = p_new
        pi = pi_new
        w = w_new

        # Update concentration parameters (prior Beta distribution)
        a, b = beta_parameters_from_beta_values(p)

        # Compute log-likelihood
        loss = np.mean(-beta_ll(p, a, b) - binomial_ll(n, k, pi))
        pbar.set_description(f'Infer beta values (NLL={loss})')

        # Stopping criterion
        losses.append(loss)
        #if abs((losses[-2] - losses[-1]) / losses[-1]) < tol:
        #    break

    return p, a, b, w


def beta_values_from_beta_parameters_and_confounder(
        n: np.ndarray,
        k: np.ndarray,
        a: np.ndarray,
        b: np.ndarray,
        y: np.ndarray,
        w: np.ndarray,
        eps: float = 1e-7
) -> np.ndarray:

    # Compute initial solution
    alpha = k + a
    beta = n - k + b
    p = alpha / (alpha + beta)
    p = np.clip(p, eps, 1 - eps)

    # Maximize full log-likelihood function
    p = solve_bayesian_logistic_regressions_wrt_p(n, k, p, y, w, a, b, eps=eps)

    return np.clip(p, eps, 1 - eps)


class BetaBinomialTransformer:

    def __init__(self, tol: float = 1e-6, max_n_iter: int = 100, gamma: float = 0.1):
        self.tol: float = tol
        self.max_n_iter: int = max_n_iter
        self.gamma: float = gamma
        self.a: Optional[np.ndarray] = None
        self.b: Optional[np.ndarray] = None
        self.w: Optional[np.ndarray] = None

    def fit(self, n: np.ndarray, k: np.ndarray, y: Optional[np.ndarray] = None) -> 'BetaBinomialTransformer':
        if y is not None:
            p, a, b, w = infer_beta_values_with_confounders(n, k, y, tol=self.tol, max_n_iter=self.max_n_iter, gamma=self.gamma)
            self.w = w
        else:
            p, a, b = infer_beta_values(n, k, tol=self.tol, max_n_iter=self.max_n_iter)
        self.a = a
        self.b = b
        return self

    def fit_transform(self, n: np.ndarray, k: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        if y is not None:
            p, a, b, w = infer_beta_values_with_confounders(n, k, y, tol=self.tol, max_n_iter=self.max_n_iter, gamma=self.gamma)
            self.w = w
        else:
            p, a, b = infer_beta_values(n, k, tol=self.tol, max_n_iter=self.max_n_iter)
        self.a = a
        self.b = b
        return p

    def transform(self, n: np.ndarray, k: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        if y is not None:
            return beta_values_from_beta_parameters_and_confounder(n, k, self.a, self.b, y, self.w)
        else:
            return beta_values_from_beta_parameters(n, k, self.a, self.b)
