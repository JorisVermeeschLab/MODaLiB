from typing import Optional

import numpy as np
import scipy.special
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import RobustScaler

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

        alpha = SemiSupervisedDeconvolution(1).fit_transform([X], [F_PROFILES])
        alpha = alpha[:, [2, 1, 0]]
        alpha /= np.sum(alpha, axis=1)[:, np.newaxis]
        return alpha

        import scipy.optimize
        alpha = np.zeros((len(X), 3))
        for i in range(len(X)):
            y, _ = scipy.optimize.nnls(F_PROFILES.T, X[i, :])
            y /= np.sum(y)
            print(y)
            alpha[i, :] = y[[2, 1, 0]]
        #alpha /= np.sum(alpha, axis=1)[:, np.newaxis]
        return alpha

        """
        # GC-content
        gc = np.zeros(len(data['Coverage']))
        for i in range(len(data['Coverage'])):
            n = data['Coverage'][i]
            gc[i] = np.sum(n * data['GC-content'][i]) / np.sum(n)
        print(gc)
        
        # Methylation
        methylation = np.zeros(len(data['CpG count']))
        for i in range(len(data['CpG count'])):
            n = data['CpG count'][i]
            methylation[i] = np.sum(n * data['CpG Methylation (first moment)'][i]) / np.sum(n)

        # GC-correction of CpG methylation
        coef = np.polyfit(gc, methylation, 1)
        methylation -= np.polyval(coef, gc)
        
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

        # Nucleosome positioning
        nuc = np.zeros(len(data['Coverage']))
        for i in range(len(data['Coverage'])):
            n = data['Coverage'][i]
            nuc[i] = np.sum(n * data['Nucleosome positioning'][i]) / np.sum(n)

        # Features
        features = np.asarray([
            methylation,
            emh['CGNN'] / emh['NCGN'],
            (emh['ANNC'] + emh['ANNT'] + emh['GNNC'] + emh['GNNT']),  #emh['ANNC'],
            (emh['CGNN'] + emh['TGNN']),
            (emh['CCNN']),
            JaggednessModality.compute_jaggedness(data),
            subnucleosomal,
            nucleosomal,
            oligonucleosomal,
            nuc
        ]).T
        if fit:
            features = self.scaler.fit_transform(features).T
        else:
            features = self.scaler.transform(features).T


        import seaborn
        import matplotlib.pyplot as plt

        P = np.zeros((features.T.shape[1], y.shape[1]))
        C = np.zeros((features.T.shape[1], y.shape[1]))
        for i in range(features.T.shape[1]):
            for j in range(y.shape[1]):
                res = scipy.stats.spearmanr(features.T[:, i], y[:, j])
                C[i, j] = res.statistic
                P[i, j] = -np.log10(res.pvalue)
        seaborn.heatmap(P, annot=C)
        plt.xticks(np.arange(3) + 0.5, ['DFFB', 'DNASE1', 'DNASE1L3'])
        plt.yticks(np.arange(features.T.shape[1]) + 0.5, ['Meth', 'CGN/NCG', 'A/G', 'CG/TG', 'CC', 'Jaggedness', 'Subnucl', 'Nucl', 'Oligonucl', 'NUC'], rotation=0)
        plt.show()

        return
        


        features[2:5, :] *= 1. / 3.
        features[6:9, :] *= 1. / 3.

        R_sign = np.asarray([
            [-1, -1,  1, -1, -1, -1, -1, -1,  1,  1], # dffb
            [-1, -1, -1,  1, -1,  1,  1, -1, -1,  1], # dnase1
            [ 1,  1, -1, -1,  1,  1, -1,  1, -1, -1]  # dnase1l3
        ])

        #features = features[:5, :]
        #R_sign = R_sign[:, :5]

        if fit:
            nuclease_activities, R = constrained_sign_deconvolution(features.T, R_sign, max_n_iter=self.max_n_iter)
            self.R = R
        else:
            nuclease_activities, _, _ = constrained_sign_deconvolution(features.T, R_sign, max_n_iter=self.max_n_iter, R=self.R)
    
        # Normalize activities
        nuclease_activities = np.log(nuclease_activities)
        nuclease_activities *= 0.4
        nuclease_activities = scipy.special.softmax(nuclease_activities, axis=1)

        return 100 * nuclease_activities
        """
