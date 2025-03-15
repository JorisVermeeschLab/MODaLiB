from typing import Optional

import numpy as np
import scipy.special
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import RobustScaler, StandardScaler

from .modalities.base import MultimodalData
from .modalities.jaggedness import JaggednessModality
from .deconvolution.constrained_sign import constrained_sign_deconvolution, unconstrained_sign_deconvolution
from .deconvolution.semi_supervised import SemiSupervisedDeconvolution
from .end_motifs import EndMotifHistograms, F_PROFILES


class NucleaseActivityEstimator(BaseEstimator, TransformerMixin):

    def __init__(self, max_n_iter: int = 5000):
        self.scaler: RobustScaler = RobustScaler()
        self.max_n_iter: int = max_n_iter
        self.R: Optional[np.ndarray] = None

    def fit(self, data: MultimodalData, y: Optional[np.ndarray] = None):
        self.run(data, y=y, fit=True)
        return self

    def fit_transform(self, data: MultimodalData, y: Optional[np.ndarray] = None) -> np.ndarray:
        return self.run(data, y=y, fit=True)

    def transform(self, data: MultimodalData) -> np.ndarray:
        return self.run(data, fit=False)

    def run(self, data: MultimodalData, y: Optional[np.ndarray] = None, fit: bool = True) -> np.ndarray:

        X = data['End-motif probability (5 prime)']
        X = X / np.sum(X, axis=1)[:, np.newaxis]
        
        # Methylation
        methylation = np.zeros(len(data['CpG count']))
        for i in range(len(data['CpG count'])):
            n = data['CpG count'][i]
            methylation[i] = np.sum(n * data['CpG Methylation (first moment)'][i]) / np.sum(n)
        
        # End motif histograms
        emh = EndMotifHistograms(data['End-motif probability (5 prime)'])

        # Fragment length distributions
        frag_length_dist = np.asarray(data['Fragment length counts'])
        frag_length_dist[:, 40:-1] /= np.sum(frag_length_dist[:, 40:-1], axis=1)[:, np.newaxis]
        subnucleosomal = np.mean(frag_length_dist[:, 40:100], axis=1)
        nucleosomal = np.mean(frag_length_dist[:, 150:250], axis=1)   #nucleosomal = np.mean(frag_length_dist[:, 150:190], axis=1)
        oligonucleosomal = np.mean(frag_length_dist[:, 300:-1], axis=1)
        norm = subnucleosomal + nucleosomal + oligonucleosomal
        subnucleosomal /= norm
        nucleosomal /= norm
        oligonucleosomal /= norm

        # Features
        features = np.asarray([
            emh['CGNN'] / emh['NCGN'],
            emh['ANNC'] + emh['ANNT'] + emh['GNNC'] + emh['GNNT'],
            emh['CGNN'] + emh['TGNN'],
            emh['CCNN'],
            JaggednessModality.compute_jaggedness(data),
            subnucleosomal,
            nucleosomal,
            oligonucleosomal
        ]).T
        if fit:
            features = self.scaler.fit_transform(features).T
        else:
            features = self.scaler.transform(features).T

        features[1:4, :] *= 1. / 3.
        features[5:8, :] *= 1. / 3.

        R_sign = np.asarray([
            [-1,  1, -1, -1, -1, -1, -1,  1], # dffb
            [-1, -1,  1, -1,  1,  1, -1, -1], # dnase1
            [ 1, -1, -1,  1,  1, -1,  1, -1]  # dnase1l3
        ])

        if fit:
            nuclease_activities, R = constrained_sign_deconvolution(features.T, R_sign, max_n_iter=self.max_n_iter)
            self.R = R
        else:
            nuclease_activities, _, _ = constrained_sign_deconvolution(features.T, R_sign, max_n_iter=self.max_n_iter, R=self.R)
    
        # Normalize activities
        nuclease_activities = np.log(nuclease_activities)
        #nuclease_activities = RobustScaler().fit_transform(nuclease_activities)
        #nuclease_activities[:, 0] -= 0.2
        #nuclease_activities[:, 2] += 0.2
        nuclease_activities *= 0.3
        nuclease_activities = scipy.special.softmax(nuclease_activities, axis=1)

        return 100 * nuclease_activities
