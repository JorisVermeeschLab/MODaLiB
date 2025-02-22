from typing import Dict, Tuple

import tqdm
import numpy as np
import scipy.stats
import pandas as pd
import statsmodels.api
from sklearn.preprocessing import OneHotEncoder, robust_scale
from statsmodels.nonparametric.smoothers_lowess import lowess


def build_design_matrix(confounders: Dict[str, Tuple[np.ndarray, str]]) -> np.ndarray:
    X = []
    for var_name in confounders.keys():
        values, var_type = confounders[var_name]
        assert var_type in {'continuous', 'binary', 'categorical'}
        assert len(values.shape) == 1
        if var_type == 'continuous':
            # If continuous variable, standardize it
            x = np.squeeze(robust_scale(values[:, np.newaxis]))[:, np.newaxis]
        elif var_type == 'binary':
            x = values[:, np.newaxis]
        else:
            x = OneHotEncoder(sparse=False).fit_transform(values[:, np.newaxis])
            x = x[:, :-1]  # Remove last column to avoid redundancy
        X.append(x)

    # Add intercept
    n = len(X[-1])
    X.append(np.ones((n, 1)))

    return np.concatenate(X, axis=1)


def bias_correction(Y_success: np.ndarray, Y_failure: np.ndarray, gc_content: np.ndarray, confounders: Dict[str, Tuple[np.ndarray, str]]) -> np.ndarray:

    assert(len(Y_success.shape) == 2)
    assert(Y_success.shape == Y_failure.shape)

    Y_corrected = np.empty_like(Y_success)

    # X = build_design_matrix(confounders)

    start_params = None
    for k in tqdm.tqdm(range(Y_success.shape[1])):  # Iterate over explained variables

        # Generalized linear model for bias correction
        # Nucleosomes: 
        # Fragment lengths: Binomial
        # Methylation: Binomial
        # Fragment counts: 
        # Mutational burden: 

        beta = Y_success[:, k] / (Y_success[:, k] + Y_failure[:, k])

        X = np.concatenate((
            gc_content[:, k, np.newaxis],
            np.ones((len(gc_content), 1))
        ), axis=1)

        family = statsmodels.api.families.Binomial()
        endog = np.asarray([Y_success[:, k], Y_failure[:, k]]).T
        model = statsmodels.api.GLM(endog, X, family=family)
        try:
            res = model.fit(start_params=start_params, maxiter=100)
        except ValueError:
            continue
        start_params = res.params

        mu = model.predict(res.params)

        #import matplotlib.pyplot as plt
        #plt.scatter(beta, gc_content[:, k])
        #plt.show()

        # y_new = scipy.stats.binom.cdf(Y_success[:, k], Y_success[:, k] + Y_failure[:, k], mu)
        y_new = Y_success[:, k] / (Y_success[:, k] + Y_failure[:, k]) - mu

        Y_corrected[:, k] = y_new

    return Y_corrected


def coverage_bias_correction(Y: np.ndarray, gc_content: np.ndarray, confounders: Dict[str, Tuple[np.ndarray, str]]) -> np.ndarray:

    assert(len(Y.shape) == 2)
    assert(Y.shape == gc_content.shape)

    Y_corrected = np.empty_like(Y)

    # X = build_design_matrix(confounders)

    start_params = None
    for k in tqdm.tqdm(range(Y.shape[1])):  # Iterate over explained variables

        X = np.concatenate((
            gc_content[:, k, np.newaxis],
            np.ones((len(gc_content), 1))
        ), axis=1)

        family = statsmodels.api.families.Gaussian()
        endog = Y[:, k]
        model = statsmodels.api.GLM(endog, X, family=family)
        try:
            res = model.fit(start_params=start_params, maxiter=100)
        except ValueError:
            continue
        start_params = res.params

        mu = model.predict(res.params)

        #import matplotlib.pyplot as plt
        #plt.scatter(beta, gc_content[:, k])
        #plt.show()

        # y_new = scipy.stats.binom.cdf(Y_success[:, k], Y_success[:, k] + Y_failure[:, k], mu)

        # Y_corrected[:, k] = Y[:, k] - mu

    return Y_corrected


def gc_correction(Y: np.ndarray, gc_content: np.ndarray) -> np.ndarray:
    Y_corrected = np.empty_like(Y)
    for i in tqdm.tqdm(range(len(Y)), desc='GC-correction'):
        Y_corrected[i, :] = Y[i, :] / lowess(Y[i, ::20], gc_content[i, ::20], xvals=gc_content[i, :], it=1, frac=0.3, return_sorted=False)
        # Y_corrected[i, :] = Y[i, :] / lowess(Y[i, :], gc_content[i, :], frac=0.3, return_sorted=False)
    return Y_corrected

