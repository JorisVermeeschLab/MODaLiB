import tqdm
import numpy as np
#from statsmodels.nonparametric.smoothers_lowess import lowess


def gc_correction(X, gc, order: int = 9) -> np.ndarray:
    X_out = np.copy(X)
    # for i in tqdm.tqdm(range(len(X_)), desc='GC-correction'):
    for i in tqdm.tqdm(range(len(X)), desc='GC-correction'):
        xs = gc[i, :]
        xs = (xs - np.mean(xs)) / np.std(xs)

        ys = X[i, :]
        coef = np.polyfit(xs, ys, order)
        ys_hat = np.polyval(coef, xs)

        X_out[i, :] = ys / ys_hat
        #X_out[i, :] = ys / lowess(ys, gc[i, :], it=1, frac=0.3, return_sorted=False)

    return X_out
