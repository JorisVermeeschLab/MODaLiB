from typing import List, Dict
import copy

import numpy as np
from sklearn.svm import SVC, SVR

from ..modalities import Modalities
from ..dataset import Dataset
from ..svm import SVM, CARSVM, MultiStageSVM
from ..stages.base import Annotation
from .mkl import KernelMatrixExtractor, MultipleKernelLearning


class Pipeline(object):

    def __init__(
            self,
            modality_enums: List[Modalities],
            model_name: str = 'svm-reimpl',
            C: float = 1.0,
            supervised_mkl: bool = False,
            random_state: int = 0xCAFE,
            platt_scaling: bool = False,
            mod_options: dict = {}
    ):
        self.modality_enums = modality_enums
        self.mod_options: dict = copy.deepcopy(mod_options)
        self.mod_options['cache'] = False
        self.extractor = KernelMatrixExtractor(modality_enums, mod_options=mod_options)
        self.mkl = MultipleKernelLearning(supervised=supervised_mkl)

        if model_name == 'stage-svm':
            self.model = MultiStageSVM(C=C, C_aux=0.1*C, random_state=random_state, platt_scaling=platt_scaling)
        elif model_name == 'svm-reimpl':
            self.model = SVM(C=C, random_state=random_state, platt_scaling=platt_scaling)
        elif model_name == 'svm':
            self.model = SVC(kernel='precomputed', probability=platt_scaling, C=C, class_weight='balanced', tol=1e-6, random_state=random_state)
        elif model_name == 'carsvm':
            self.model = CARSVM(C_c=C, C_r=C*0.01, epsilon=0.4, platt_scaling=platt_scaling)
        else:
            raise NotImplementedError(f'Unknown model "{model_name}"')
        self.model_name: str = model_name

    def fit(self, dataset: Dataset) -> 'Pipeline':

        annotations = dataset.annotations
        groups = dataset.groups
        y = dataset.y

        tf1 = dataset.tf / np.std(dataset.tf)
        ichor = np.asarray(dataset.ichor, dtype=float)
        tf2 = ichor / np.std(ichor)
        tf = 0.5 * (tf1 + tf2)

        self.extractor.preprocess(dataset.data)
        Ks = self.extractor.fit_transform(dataset.data, y)

        K = self.mkl.fit_transform(Ks, y=y)

        weights = np.full(len(K), 1. / len(K))
        weights /= np.sum(weights)
        weights = weights[:, np.newaxis, np.newaxis]
        K = np.sum(weights * K, axis=0)

        if self.model_name == 'stage-svm':
            self.model.fit(K, annotations)
        elif self.model_name == 'svm-reimpl':
            self.model.fit(K, y)
        elif self.model_name == 'svm':
            self.model.fit(K, y)
        elif self.model_name == 'carsvm':
            self.model.fit(K, np.asarray([y, tf]).T)
        else:
            raise NotImplementedError(f'Unknown model "{self.model_name}"')
        return self

    def predict_proba(self, dataset: Dataset) -> np.ndarray:
        self.extractor.preprocess(dataset.data)
        Ks = self.extractor.transform(dataset.data)
        K = self.mkl.transform(Ks)
        return self.model.predict_proba(K)

    def predict(self, dataset: Dataset) -> np.ndarray:
        y_proba = self.predict_proba(dataset)
        return np.argmax(y_proba, axis=1)
