from typing import Tuple, List, Optional, Union, Dict, Any

import numpy as np
import scipy.optimize
from sklearn.preprocessing import RobustScaler
import torch
import tqdm

from modalibpy.deconvolution.optimizer import Optimizer
from modalibpy.manifolds.base import MatrixManifold


class SemiSupervisedDeconvolution:

    def __init__(self, n_modalities: int, reg: float = 1.0):
        self.n_modalities: int = int(n_modalities)
        self.reg: float = float(reg)
        self.Rs: List[Optional[np.ndarray]] = [None for _ in range(self.n_modalities)]
        self.manifolds: List[Optional[MatrixManifold]] = [None for _ in range(self.n_modalities)]

    @staticmethod
    def to_tensors(Rs: List[Union[np.ndarray, MatrixManifold]]) -> Tuple[List[torch.Tensor], int]:
        # Convert arrays to tensors, and get the number of samples
        n = 0
        Rs_tensors = []
        for R in Rs:
            if isinstance(R, MatrixManifold):
                if (n != 0) and (n != R.n):
                    raise ValueError(f'Atlas size inconsistency: {n} != {len(R)}')
                Rs_tensors.append(R.initialize())
            else:
                if (n != 0) and (n != len(R)):
                    raise ValueError(f'Atlas size inconsistency: {n} != {len(R)}')
                n = len(R)
                Rs_tensors.append(torch.FloatTensor(R))
        return Rs_tensors, n

    @staticmethod
    def nnls(X: np.ndarray, R: np.ndarray) -> np.ndarray:
        alpha = []
        for i in range(len(X)):
            alpha_i, _ = scipy.optimize.nnls(R.T, X[i, :])
            sum_i = np.sum(alpha_i)
            if sum_i > 0:
                alpha_i /= sum_i
            alpha.append(alpha_i)
        return np.asarray(alpha)

    @staticmethod
    def criterion(X: torch.Tensor, X_hat: torch.Tensor) -> torch.Tensor:
        dim = 0  # TODO
        norm = 2

        if norm == 2:
            ss_res = torch.mean(torch.square(X - X_hat))
            ss_tot = torch.mean(torch.square(X - torch.mean(X, dim=dim).unsqueeze(dim)))
        else:
            ss_res = torch.mean(torch.abs(X - X_hat))
            ss_tot = torch.mean(torch.abs(X - torch.mean(X, dim=dim).unsqueeze(dim)))

        return ss_res / ss_tot

    @staticmethod
    def unimodal_deconvolution(
            X: torch.Tensor,
            R: torch.Tensor,
            manifold: MatrixManifold,
            alpha: torch.Tensor,
            max_n_iter: int = 1000
    ) -> np.ndarray:
        R = torch.nn.Parameter(torch.clone(R))

        # Initialize optimizer
        optimizer = Optimizer()
        optimizer.add([R], 1e-2)

        pbar = tqdm.tqdm(range(max_n_iter))
        for _ in pbar:
            optimizer.zero_grad()
            R_ret = manifold.retraction(R)
            X_hat = torch.mm(alpha, R_ret)
            loss = SemiSupervisedDeconvolution.criterion(X, X_hat)
            loss.backward()
            optimizer.step(loss.item())
            pbar.set_description(str(loss.item()))

        return R.data

    def run(
            self,
            Xs: List[np.ndarray],
            Rs: List[Union[np.ndarray, MatrixManifold]],
            fit: bool=True,
            max_n_iter: int = 1000
    ) -> np.ndarray:
        
        assert len(Xs) == len(Rs)

        # Get manifolds
        if fit:
            manifolds = [R if isinstance(R, MatrixManifold) else None for R in Rs]
        else:
            manifolds = self.manifolds

        # Convert NumPy arrays to PyTorch tensors
        Xs, n = SemiSupervisedDeconvolution.to_tensors(Xs)
        Rs, m = SemiSupervisedDeconvolution.to_tensors(Rs)

        # Initialize proportions:
        # Supervised deconvolution using only the modalities for which
        # a reference atlas is available
        if fit:
            mask = [manifold is None for manifold in manifolds]
            if not any(mask):
                mask = [True] * len(mask)
            X = np.concatenate([Xs[i].cpu().data.numpy() for i in range(len(Xs)) if mask[i]], axis=1)
            R = np.concatenate([Rs[i].cpu().data.numpy() for i in range(len(Rs)) if mask[i]], axis=1)
        else:
            X = np.concatenate([Xs[i].cpu().data.numpy() for i in range(len(Xs))], axis=1)
            R = np.concatenate(Rs, axis=1)
        alpha0 = torch.FloatTensor(SemiSupervisedDeconvolution.nnls(X, R))
        alpha_log = torch.nn.Parameter(torch.log(alpha0 + 1e-10))

        # Initialize reference entities
        if fit:
            for i in range(len(Rs)):
                if manifolds[i] is not None:
                    #Rs[i] = SemiSupervisedDeconvolution.unimodal_deconvolution(
                    #    Xs[i], Rs[i], manifolds[i], alpha0, max_n_iter=max_n_iter
                    #)
                    pass  # TODO

        # Gather parameters
        R_parameters = []
        for i in range(len(Rs)):
            if fit and (manifolds[i] is not None):
                Rs[i] = torch.nn.Parameter(Rs[i])
                R_parameters.append(Rs[i])

        # Initialize optimizer
        optimizer = Optimizer()
        if len(R_parameters) > 0:
            optimizer.add(R_parameters, 1e-2)
        optimizer.add(alpha_log, 1e-1)

        pbar = tqdm.tqdm(range(max_n_iter))
        for _ in pbar:
            optimizer.zero_grad()

            alpha = torch.softmax(alpha_log, dim=1)

            # Deconvolution error
            losses = []
            for i, (X, R, manifold) in enumerate(zip(Xs, Rs, manifolds)):

                if manifold is not None:
                    R = manifold.retraction(R)

                X_hat = torch.mm(alpha, R)

                # Loss function
                loss = SemiSupervisedDeconvolution.criterion(X, X_hat)
                losses.append(loss)

            # Regularization
            if fit and (self.reg > 0):
                reg = torch.mean(torch.abs(torch.mean(alpha0, dim=0) - torch.mean(alpha, dim=0))) / torch.mean(torch.abs(torch.mean(alpha0, dim=0).unsqueeze(0) - alpha0))
                losses.append(self.reg * reg)

            #print(losses)

            loss = sum(losses) / len(losses)

            loss.backward()

            optimizer.step(loss.item())

            pbar.set_description(str(loss.item()))

        if fit:
            self.Rs = [R.cpu().data.numpy() for R in Rs]
            self.manifolds = manifolds

        alpha = torch.softmax(alpha_log, dim=1).cpu().data.numpy()
        return alpha

    def fit(self, Xs: List[np.ndarray], Rs: List[Optional[np.ndarray]], **kwargs):
        self.run(Xs, Rs, fit=True, **kwargs)
        return self

    def fit_transform(self, Xs: List[np.ndarray], Rs: List[Optional[np.ndarray]], **kwargs) -> np.ndarray:
        return self.run(Xs, Rs, fit=True, **kwargs)

    def transform(self, Xs: List[np.ndarray], **kwargs) -> np.ndarray:
        return self.run(Xs, self.Rs, fit=False, **kwargs)
