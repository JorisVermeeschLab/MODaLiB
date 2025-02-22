import numpy as np
import scipy.spatial
import scipy.stats


class BaseStep:

    def compare(self, X, dist, dist_transformed):

        D2 = scipy.spatial.distance.cdist(X, X, metric='sqeuclidean')

        D2_hat_old = dist.sqeuclidean_cdist()
        D2_hat_new = dist_transformed.sqeuclidean_cdist()

        idx = np.triu_indices(len(D2), k=1)
        d2_hat_old = D2_hat_old[idx[0], idx[1]]
        d2_hat_new = D2_hat_new[idx[0], idx[1]]
        d2 = D2[idx[0], idx[1]]

        max_ = np.maximum(d2, d2_hat_old)
        explained_ratio_before = np.mean((max_ - np.abs(d2 - d2_hat_old)) / max_)
        max_ = np.maximum(d2, d2_hat_new)
        explained_ratio_after = np.mean((max_ - np.abs(d2 - d2_hat_new)) / max_)
        print('EXPLRATIOS', explained_ratio_before, explained_ratio_after)

        pearson_before = scipy.stats.pearsonr(d2, d2_hat_old)
        pearson_after = scipy.stats.pearsonr(d2, d2_hat_new)
        print('Pearson', pearson_before, pearson_after)

        return True  # explained_ratio_after > explained_ratio_before
