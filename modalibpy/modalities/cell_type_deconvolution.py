from typing import Optional

import numpy as np
from sklearn.preprocessing import RobustScaler
import scipy.optimize

from .base import BaseModality, MultimodalData
from ..preprocessing import BetaBinomialTransformer
from ..decorrelate import gc_correction
from ..deconvolution.semi_supervised import SemiSupervisedDeconvolution
from ..deconvolution.gene_expression import crossdomain_deconvolution
from ..pipeline.fs import FeatureSelector
from ..manifolds import HistogramMatrixManifold, RatioMatrixManifold
from ..end_motifs import cgnngc_ratios
from ..deconvolution.metdecode import MetDecode


class FragCellTypeDeconvolutionModality(BaseModality):


    def __init__(self, *args, max_n_iter: int = 2000, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='chi2', **kwargs)
        self.max_n_iter: int = max_n_iter
        self.use_compartments: bool = kwargs.get('use_compartments', False)

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('Frag cell type proportions' in data) and self.cache:
            return None

        X = data['Fragmentation-based-gene-expression']
        R = data['Gene-expression-reference']

        X = crossdomain_deconvolution(X, R, max_n_iter=self.max_n_iter)

        if self.use_compartments:
            alpha = np.zeros((len(X), 4))
            for j, compartment in enumerate(['immune', 'endothelial', 'epithelial', 'stromal']):
                mask = (np.asarray(data['Gene-expression-reference-compartments'], dtype=object) == compartment)
                alpha[:, j] = np.sum(X[:, mask], axis=1)
            X = alpha

        self.save(data, 'Frag cell type proportions', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        X = self.load(data, 'Frag cell type proportions')

        X = X + 1e-7

        X = X / np.sum(X, axis=1)[:, np.newaxis]

        return X


class MethCellTypeDeconvolutionModality(BaseModality):

    def __init__(self, *args, max_n_iter: int = 10, **kwargs) -> None:
        BaseModality.__init__(self, *args, kernel='chi2', **kwargs)
        self.scaler: RobustScaler = RobustScaler()
        self.bb_transformer: BetaBinomialTransformer = BetaBinomialTransformer(gamma=0.01, lr=1e-6)
        self.max_n_iter: int = max_n_iter
        self.gamma_hat: Optional[np.ndarray] = None

    def preprocess(self, data: MultimodalData, fit: bool = True) -> None:

        if ('Cell type proportions' in data) and self.cache:
            return None

        X = data['CpG Methylation (deconvolution, cfDNA)']
        R = data['CpG Methylation (deconvolution, reference)']

        y_cgn_ncg = cgnngc_ratios(data['End-motif probability (5 prime)'])
        Y = np.repeat(y_cgn_ncg[:, np.newaxis], X.shape[1], axis=1)

        n = np.clip(data['CpG count (deconvolution, cfDNA)'], 0.01, None)

        n_atlas_entities = len(data['CpG Methylation (deconvolution, reference)'])

        if self.probabilistic:
            if fit:
                X = self.bb_transformer.fit_transform(n, X * n, y=Y)
            else:
                X = self.bb_transformer.transform(n, X * n, y=Y)

        # Supervised deconvolution
        use_coverage = False
        n_unknowns = 0
        unsupervised = False
        if fit:
            model = MetDecode(
                R * 10000,
                np.ones_like(R) * 10000,  # Atlas is assumed to have sufficient coverage
                X * n,
                n,
                n_unknowns,
                coverage=use_coverage,
                unsupervised=unsupervised,
                sum1=True
            )
            X = model.deconvolute()
            self.gamma_hat = model.gamma_hat
        else:
            model = MetDecode(
                self.gamma_hat * 10000,
                np.ones_like(R) * 10000,  # Atlas is assumed to have sufficient coverage
                X * n,
                n,
                n_unknowns,
                coverage=use_coverage,
                unsupervised=False,
                sum1=True
            )
            X = model.deconvolute()

        self.save(data, 'Cell type proportions', X)

    def process(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:
        X = self.load(data, 'Cell type proportions')

        X = X / np.sum(X, axis=1)[:, np.newaxis]

        return X
