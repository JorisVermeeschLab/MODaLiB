import tqdm
import numpy as np
import scipy.spatial

from modalibpy.pipeline.base import BaseStep


class Decorrelate(BaseStep):

    def __init__(self):
        BaseStep.__init__(self)

    def fit_transform(self, dist, endog, exog, order: int = 9) -> np.ndarray:
        X_out = np.copy(endog)
        # for i in tqdm.tqdm(range(len(endog)), desc='Decorrelate'):
        for i in range(len(endog)):
            xs = exog[i, :]
            xs = (xs - np.mean(xs)) / np.std(xs)

            ys = endog[i, :]
            coef = np.polyfit(xs, ys, order)
            ys_hat = np.polyval(coef, xs)
            #X_out[i, :] = np.clip(ys / ys_hat, 0.5, 2)
            X_out[i, :] = ys_hat

        endog = endog - X_out

        dist_transformed = dist - X_out

        if not self.compare(endog, dist, dist_transformed):
            dist_transformed = dist

        return dist_transformed, endog
