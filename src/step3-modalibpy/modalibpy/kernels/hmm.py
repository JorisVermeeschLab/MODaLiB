from typing import Optional, List

import tqdm
import numpy as np
import scipy.spatial

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import RobustScaler

from .base import Kernel


class HMMKernel(BaseEstimator, TransformerMixin):

    def __init__(self):
        from hmmlearn.hmm import GaussianHMM
        from hmmlearn import _hmmc
        self.X_train: Optional[np.ndarray] = None
        self.hmm: Optional[GaussianHMM] = None
        self.gamma: float = 1

    @staticmethod
    def train_hmm(X: np.ndarray):
        n_states = 5
        e = 0.99999
        I = np.eye(n_states)
        transmat_prior = I * e + (1 - I) * ((1 - e) / (n_states - 1))
        transmat_prior *= 100

        n = X.shape[0]
        p = X.shape[1]

        hmm = GaussianHMM(n_components=n_states, transmat_prior=transmat_prior, init_params='mcs')
        hmm.fit(X.flatten()[:, np.newaxis], lengths=np.full(n, p, dtype=int))
        #hmm.fit(X[i, :, np.newaxis])
        #if np.any(np.sum(hmm.transmat_, axis=1) == 0):
        #    hmm = GaussianHMM(n_components=n_states, init_params='mcs')
        #    hmm.fit(X[i, :, np.newaxis])

        # Sort states by mean
        #means = np.squeeze(hmm.means_)
        #idx = np.argsort(means)
        #hmm.means_ = hmm.means_[idx, :]
        #hmm.covars_ = hmm.covars_[idx, :].reshape(n_states, 1)
        #hmm.startprob_ = hmm.startprob_[idx]
        #hmm.transmat_ = hmm.transmat_[idx, :][:, idx]

        return hmm

    @staticmethod
    def _compute_fisher_gradient(hmm, x: np.ndarray) -> np.ndarray:

        print('k5')

        means = np.squeeze(hmm.means_)[np.newaxis, :]
        variances = np.squeeze(hmm.covars_)[np.newaxis, :]

        n_obs = len(x)
        n_states = len(means)

        eps = 1e-7

        startprob = np.clip(hmm.startprob_, eps, None)
        startprob = startprob / np.sum(startprob)

        transmat = np.clip(hmm.transmat_, eps, None)  # best: 1e-4
        transmat = transmat / np.sum(transmat, axis=1)[:, np.newaxis]

        gamma = hmm.predict_proba(x[:, np.newaxis])

        mean_gradient = -2 * (x[:, np.newaxis] - means) / variances
        #mean_gradient = np.mean(gamma * mean_gradient, axis=0)
        mean_gradient = np.sum(gamma * mean_gradient, axis=0)

        #variance_gradient = np.mean(gamma * variance_gradient, axis=0)
        variance_gradient = -0.5 * np.sum(gamma / variances, axis=0)
        variance_gradient += 0.5 * np.sum(gamma * np.square((x[:, np.newaxis] - means) / variances), axis=0)
        #variance_gradient = 0.5 * np.sum(gamma * np.square((x[:, np.newaxis] - means) / variances), axis=0)

        # Gradient w.r.t. start probabilities
        gamma0 = np.clip(gamma[0, :], eps, None)
        start_prob_gradient = gamma0 / startprob

        # Gradient w.r.t. transition probabilities
        frameprob = hmm._compute_likelihood(x[:, np.newaxis])
        log_prob, fwdlattice, scaling_factors = _hmmc.forward_scaling(hmm.startprob_, transmat, frameprob)
        bwdlattice = _hmmc.backward_scaling(hmm.startprob_, transmat, frameprob, scaling_factors)
        xi_sum = _hmmc.compute_scaling_xi_sum(fwdlattice, transmat, bwdlattice, frameprob)
        xi_sum = np.clip(xi_sum, eps, None)
        transmat_gradient = (xi_sum / transmat).flatten()

        gradient = np.concatenate((mean_gradient, variance_gradient, start_prob_gradient, transmat_gradient), axis=0)
        #gradient = np.concatenate((startprob, transmat.flatten(), means.flatten(), variances.flatten()), axis=0)

        #gradient *= 0.0001

        #gradient = gradient / np.linalg.norm(gradient)

        return gradient

    def compute_fisher_gradients(self, X: np.ndarray) -> np.ndarray:

        gradients = []
        for i in tqdm.tqdm(range(len(X)), desc='Computing Fisher gradients'):
            gradients.append(HMMKernel._compute_fisher_gradient(self.hmm, X[i, :]))
        gradients = np.asarray(gradients)

        fim = np.einsum('ik,il->kl', gradients, gradients) / len(gradients)
        M = np.linalg.inv(scipy.linalg.sqrtm(fim))
        M = np.real(M)
        gradients = np.dot(gradients, M)

        #gradients = gradients / np.linalg.norm(gradients, axis=1)[:, np.newaxis]

        #gradients = RobustScaler().fit_transform(gradients)

        return gradients

    def pairwise_kernel(self, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        if X is Y:
            x_gradients = self.compute_fisher_gradients(X)
            y_gradients = x_gradients
            D2 = scipy.spatial.distance.cdist(x_gradients, x_gradients, metric='sqeuclidean')
            #D2 = scipy.spatial.distance.cdist(X, X, metric='sqeuclidean')
            self.gamma = Kernel.infer_bandwidth(D2)
        else:
            x_gradients = self.compute_fisher_gradients(X)
            y_gradients = self.compute_fisher_gradients(Y)
            D2 = scipy.spatial.distance.cdist(x_gradients, y_gradients, metric='sqeuclidean')
            #D2 = scipy.spatial.distance.cdist(X, Y, metric='sqeuclidean')
        return Kernel.safe_exp(-self.gamma * D2)
        #return np.dot(x_gradients, y_gradients.T)

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> 'RBFKernel':
        self.X_train = X
        self.hmm = HMMKernel.train_hmm(X)
        return self

    def transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        return self.pairwise_kernel(X, self.X_train)

    def fit_transform(self, X: np.ndarray, y: Optional[np.ndarray] = None) -> np.ndarray:
        self.X_train = X
        self.hmm = HMMKernel.train_hmm(X)
        K = self.pairwise_kernel(X, X)

        import matplotlib.pyplot as plt
        plt.imshow(K)
        plt.colorbar()
        plt.show()

        return K
