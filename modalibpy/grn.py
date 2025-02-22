from typing import List, Dict, Tuple

import numpy as np
import scipy.stats
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score


class GRN(object):

    def __init__(self, M: np.ndarray, tfs: List[str], tgs: List[str]):
        self.M: np.ndarray = M
        self.tfs: List[str] = tfs
        self.tgs: List[str] = tgs

    def to_numpy(self) -> np.ndarray:
        return self.M

    @staticmethod
    def load(filepath: str) -> 'GRN':
        df = pd.read_csv(filepath, sep=',', header='infer')
        tfs = list(set(df['source'].to_numpy()))
        tgs = list(set(df['target'].to_numpy()))
        tf_dict = {name: i for i, name in enumerate(tfs)}
        tg_dict = {name: i for i, name in enumerate(tgs)}
        M = np.zeros((len(tfs), len(tgs)), dtype=float)
        for tf, tg, weight in zip(df['source'], df['target'], df['weight']):
            i = tf_dict[tf]
            j = tg_dict[tg]
            M[i, j] = weight
        return GRN(M, tfs, tgs)

    def reorder(self, tfs: List[str], tgs: List[str]) -> 'GRN':
        tf_dict = {name: i for i, name in enumerate(self.tfs)}
        tg_dict = {name: i for i, name in enumerate(self.tgs)}
        idx_tf = np.asarray([tf_dict[tf] for tf in tfs], dtype=int)
        idx_tg = np.asarray([tg_dict[tg] for tg in tgs], dtype=int)
        M = self.M[idx_tf, :][:, idx_tg]
        return GRN(M, tfs, tgs)

    def score(self, grn_pred: 'GRN') -> Dict[str, float]:
        shared_tfs = list(set(self.tfs).intersection(grn_pred.tfs))
        shared_tgs = list(set(self.tgs).intersection(grn_pred.tgs))
        goldstandard = self.reorder(shared_tfs, shared_tgs)
        grn_pred = grn_pred.reorder(shared_tfs, shared_tgs)
        y_true = (np.abs(goldstandard.to_numpy().flatten()) > 0).astype(int)
        y_score = np.abs(grn_pred.to_numpy().flatten())
        print(np.mean(y_true))
        print(y_score)
        return {
            'auroc': roc_auc_score(y_true, y_score),
            'auprc': average_precision_score(y_true, y_score),
            'n-tfs': len(shared_tfs),
            'n-tgs': len(shared_tgs),
            't-test': scipy.stats.ttest_ind(y_score[y_true == 1], y_score[y_true == 0]).pvalue
        }
