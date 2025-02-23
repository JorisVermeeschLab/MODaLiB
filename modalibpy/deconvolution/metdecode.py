from typing import Tuple, Any, Optional

import tqdm
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.special
from scipy.optimize import nnls
from scipy.special import logit
import torch
import torch.nn.functional

from .optimizer import Optimizer


def chr_name_to_int(chr_name: str) -> int:
    """Convert chromosome name to identifier.

    Args:
        chr_name: Chromosome name (e.g., "chr8").

    Returns:
        Integer representation of the chromosome.
    """
    chr_name = chr_name.lower()
    if chr_name.startswith('chr'):
        chr_name = chr_name[3:]
        if chr_name == 'x':
            chr_id = 23
        elif chr_name == 'y':
            chr_id = 24
        else:
            chr_id = int(chr_name)
    else:
        chr_id = int(chr_name)
    return chr_id - 1


def row_name_as_int(combined_coords: str) -> int:
    """Convert genomic coordinates to integer.

    These integers are further used to sort marker regions
    by genomic coordinates efficiently.

    Args:
        combined_coords: Concatenation of genomic coordinates
            of the form "chr1-3792834-3793456".

    Returns:
        Integer representation of the marker region.
    """
    elements = combined_coords.split('-')
    assert len(elements) == 3
    chr_name = elements[0]
    start = int(elements[1])
    return 1000000000 * chr_name_to_int(chr_name) + start


def load_input_file(
        filepath: str,
        chr_idx: Optional[np.ndarray] = None,
        values_idx: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load input TSV file containing the counts.

    Each tissue / cell type has two dedicated columns, namely the number of methylated CpG
    sites spanned in the marker region, and the total number of CpG sites
    (both methylation and unmethylated). Each row corresponds to a marker region.
    The first 3 columns contain respectively the chromosome, start position and end position
    of each marker region. The file must contain a header of the form:
    CHROM    START   END TISSUE1_METH    TISSUE1_DEPTH   TISSUE2_METH    ...

    Args:
        filepath: Input TSV file.

    Returns:
        methylated: A matrix of shape `(n_samples, n_markers)` containing the methylated counts
            for each sample and marker region.
        depths: A matrix of shape `(n_samples, n_markers)` containing the total counts
            (methylated + unmethylated) for each sample and marker region.
        column_names: List of column (marker) names.
        row_names: List of row (sample) names.
    """

    # Check the presence of the header
    with open(filepath, 'r') as f:
        line = f.readline()
        if line.startswith('CHROM'):
            header = 'infer'
        else:
            header = None

    # Parse input file
    df = pd.read_csv(filepath, delimiter='\t', header=header)

    # Reindex profiles by genomic coordinates
    if chr_idx is None:
        row_names = df.iloc[:, 0].astype(str) + '-' + df.iloc[:, 1].astype(str) + '-' + df.iloc[:, 2].astype(str)
    else:
        row_names = ''
        for k, i in enumerate(chr_idx):
            row_names = row_names + df.iloc[:, i].astype(str)
            if k < len(chr_idx) - 1:
                row_names = row_names + '-'
    df = df.set_index(row_names, drop=False)

    # Retrieve the count matrices
    if values_idx is None:
        values = df.iloc[:, 3:].values.T
    else:
        values = df.iloc[:, values_idx].values.T
    df.reset_index()
    methylated = values[::2, :].astype(float)
    depths = values[1::2, :].astype(int)

    # Row and column labels
    if header is not None:
        if values_idx is None:
            column_names = [s.replace('_METH', '') for s in df.columns[3::2]]
        else:
            column_names = [s.replace('_METH', '') for s in df.columns[values_idx][::2]]
    else:
        column_names = [f'sample-{j + 1}' for j in range(depths.shape[0])]
    column_names = np.asarray(column_names, dtype=object)
    row_names = np.asarray(row_names, dtype=object)

    # Sort marker regions by increasing genomic coordinates
    idx = np.argsort(np.asarray([row_name_as_int(row_name) for row_name in row_names], dtype=np.uint64))
    row_names = row_names[idx]
    methylated = methylated[:, idx]
    depths = depths[:, idx]

    # Ensure methylation ratios are within bounds [0, 1]
    methylated = np.clip(methylated, 0, depths)

    return methylated, depths, column_names, row_names


class MetDecode:

    def __init__(
            self,
            R_methylated: np.ndarray,
            R_depths: np.ndarray,
            X_methylated: np.ndarray,
            X_depths: np.ndarray,
            n_unknown_tissues: int = 0,
            coverage: bool = True,
            unsupervised: bool = True,
            sum1: bool = False
    ):

        self.X_methylated: np.ndarray = X_methylated
        self.X_depths: np.ndarray = X_depths
        self.R_methylated: np.ndarray = R_methylated
        self.R_depths: np.ndarray = R_depths
        self.n_known_tissues: int = self.R_methylated.shape[0]
        self.n_unknown_tissues: int = n_unknown_tissues
        self.n_tissues: int = self.n_known_tissues + self.n_unknown_tissues
        self.coverage: bool = coverage
        self.unsupervised: bool = unsupervised
        self.sum1: bool = sum1

        # Add pseudo-counts
        self.R_methylated, self.R_depths = MetDecode.add_pseudo_counts(self.R_methylated, self.R_depths)
        self.X_methylated, self.X_depths = MetDecode.add_pseudo_counts(self.X_methylated, self.X_depths)

        R = np.asarray([self.R_methylated, self.R_depths - self.R_methylated])
        self.R: np.ndarray = np.transpose(R, (1, 2, 0)).reshape(self.R_depths.shape[0], 2 * self.R_depths.shape[1])
        X = np.asarray([self.X_methylated, self.X_depths - self.X_methylated])
        self.X: np.ndarray = np.transpose(X, (1, 2, 0)).reshape(self.X_depths.shape[0], 2 * self.X_depths.shape[1])

        self.info = {
            'loss': []
        }

    @staticmethod
    def add_pseudo_counts(methylated, depths, pc=0.8):
        methylated = np.asarray(methylated, dtype=float)
        depths = np.asarray(depths, dtype=float)
        avg_meth = np.sum(methylated) / np.sum(depths)

        # Compute per-row methylation
        row_methylated = np.sum(methylated, axis=1)
        row_depths = np.sum(depths, axis=1)
        mask = np.logical_or(row_methylated == 0, row_methylated == row_depths)
        row_depths[mask] += pc
        row_methylated[mask] += pc * avg_meth
        row_meth = row_methylated / row_depths

        # Compute per-column methylation
        col_methylated = np.sum(methylated, axis=0)
        col_depths = np.sum(depths, axis=0)
        mask = np.logical_or(col_methylated == 0, col_methylated == col_depths)
        col_depths[mask] += pc
        col_methylated[mask] += pc * avg_meth
        col_meth = col_methylated / col_depths

        # Impute missing values
        meth_prior = np.outer(np.ones(len(row_meth)), col_meth)  # TODO: use survival function of Beta distribution instead?
        # TODO: would be better taking the sqrt
        meth_prior = np.clip(meth_prior, 0.01, 0.99)

        mask = np.logical_or(methylated == 0, methylated == depths)
        methylated[mask] += pc * meth_prior[mask]
        depths[mask] += pc

        assert np.all(methylated > 0)
        assert not np.any(np.isnan(methylated))
        assert not np.any(np.isnan(depths))
        assert not np.any(np.isinf(methylated))
        assert not np.any(np.isinf(depths))
        assert np.all(methylated < depths)

        return methylated, depths

    @staticmethod
    def nnls(R_atlas: np.ndarray, R_cfdna: np.ndarray, W_cfdna: np.ndarray) -> np.ndarray:
        n_samples = len(R_cfdna)
        alpha_hat = []
        for i in range(n_samples):
            A = (np.sqrt(W_cfdna[np.newaxis, i, :]) * R_atlas).T
            b = np.sqrt(W_cfdna[i, :]) * R_cfdna[i, :]
            # x, residuals = nnls(A, b, maxiter=3000)
            x, residuals = nnls(A, b)
            alpha_hat.append(x)
        alpha_hat = np.asarray(alpha_hat)
        alpha_hat /= alpha_hat.sum(axis=1)[:, np.newaxis]
        return alpha_hat

    def _init_alpha_and_gamma(self) -> Tuple[np.ndarray, np.ndarray]:
        R_atlas = self.R_methylated / self.R_depths
        R_cfdna = self.X_methylated / self.X_depths
        W_cfdna = self.X_depths / np.sum(self.X_depths)
        alpha_hat = MetDecode.nnls(R_atlas, R_cfdna, W_cfdna)
        alpha_hat += 1e-6
        alpha_hat /= np.sum(alpha_hat, axis=1)[:, np.newaxis]

        #lb = np.min(R_atlas, axis=0)
        #ub = np.max(R_atlas, axis=0)
        lb = np.quantile(R_atlas, 0.4, axis=0)
        ub = np.quantile(R_atlas, 0.6, axis=0)

        n_unknowns = self.n_unknown_tissues
        while n_unknowns > 0:
            residuals = np.dot(alpha_hat, R_atlas) - R_cfdna
            residuals = np.median(residuals, axis=0)
            r = (residuals <= 0).astype(float)
            r = lb + r * (ub - lb)
            R_atlas = np.concatenate((R_atlas, r[np.newaxis, :]), axis=0)

            alpha_hat = MetDecode.nnls(R_atlas, R_cfdna, W_cfdna)
            alpha_hat += 1e-6
            if self.sum1:
                alpha_hat /= np.sum(alpha_hat, axis=1)[:, np.newaxis]

            n_unknowns -= 1

        alpha_logit = np.log(alpha_hat)
        gamma_logit = logit(R_atlas)

        return alpha_logit, gamma_logit

    def deconvolute(self, max_n_iter: int = 2000, patience: int = 1000) -> np.ndarray:

        self.info['loss'] = []

        _alpha_logit, _gamma_logit = self._init_alpha_and_gamma()
        alpha_logit = torch.nn.Parameter(torch.FloatTensor(np.copy(_alpha_logit)))
        gamma_logit = torch.nn.Parameter(torch.FloatTensor(np.copy(_gamma_logit)))
        gamma = torch.sigmoid(gamma_logit)

        R_atlas = scipy.special.expit(_gamma_logit)
        R_cfdna = self.X_methylated / self.X_depths

        R_cfdna = torch.FloatTensor(R_cfdna)

        optimizer = Optimizer(verbose=False)
        optimizer.add([alpha_logit], 1e-2)
        if self.unsupervised:
            optimizer.add([gamma_logit], 1e-5)

        n_steps_without_improvement = 0
        best_loss = np.inf
        extra_info = {'loss': []}
        for iteration in tqdm.tqdm(range(max_n_iter)):

            optimizer.zero_grad()

            if self.sum1:
                alpha = torch.softmax(alpha_logit, dim=1)
            else:
                alpha = torch.exp(alpha_logit)
            gamma = torch.sigmoid(gamma_logit)
            X_reconstructed = torch.mm(alpha, gamma)
            if self.coverage:
                weights = torch.sqrt(torch.FloatTensor(self.X_depths))
                weights = weights / torch.mean(weights)
                loss = torch.mean(weights * torch.abs(X_reconstructed - R_cfdna))
            else:
                loss = torch.mean(torch.abs(X_reconstructed - R_cfdna))
            loss.backward()

            self.info['loss'].append(loss.item())
            extra_info['loss'].append(loss.item())

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

        alpha = scipy.special.softmax(alpha_logit.data.numpy(), axis=1)
        gamma = gamma.cpu().data.numpy()

        C = scipy.stats.spearmanr(gamma, R_atlas, axis=1)[0][:len(gamma), :][:, len(gamma):]

        A = np.zeros(C.shape)
        idx = np.argsort(np.mean(C, axis=1))
        for k1 in idx:
            k2 = np.argmax(C[k1, :])
            A[k1, k2] = 1
            C[k1, :] = -np.inf
            C[:, k2] = -np.inf

        alpha = np.dot(alpha, A)

        self.gamma_hat = np.dot(A.T, gamma)

        return alpha
