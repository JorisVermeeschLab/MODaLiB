from typing import Any, Dict, List, Union, Tuple, Optional
import traceback

import numpy as np
import scipy.stats
import tqdm
import random
from sklearn.metrics import *
from sklearn.model_selection import LeaveOneGroupOut, LeaveOneOut, KFold, StratifiedGroupKFold, GroupKFold
from sklearn.calibration import CalibratedClassifierCV
from sklearn.svm import SVC, SVR
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import KernelPCA

from .pipeline.mkl import KernelMatrixExtractor, MultipleKernelLearning
from .kernels.rbf import RBFKernel
from .modalities import Modalities
from .dataset import Dataset
from .svm import SVM, CARSVM, MultiStageSVM
from .stages.base import Annotation


def pad_predicted_proba(y_proba: np.ndarray, y_train: np.ndarray, n_classes: int) -> np.ndarray:
    n = len(y_proba)
    idx = np.asarray(list(sorted(set(y_train))), dtype=int)
    y_padded = np.zeros((n, n_classes))
    y_padded[:, idx] = y_proba
    return y_padded


def heteroscedasticity_score(K: np.ndarray, n_fragments: np.ndarray) -> float:
    pca = KernelPCA(n_components=len(K), kernel='precomputed')
    coords = pca.fit_transform(K)
    p_values = np.asarray([scipy.stats.spearmanr(coords[:, k], n_fragments).pvalue for k in range(coords.shape[1])])
    p_values = np.nan_to_num(p_values, nan=1)
    #weights = np.square(pca.eigenvalues_)
    weights = pca.eigenvalues_
    return float(np.sum(-np.log10(p_values) * weights) / np.sum(weights))


def evaluate(
        dataset: Dataset,
        modality_enums: List[Modalities],
        y: Optional[np.ndarray] = None,
        multi_stage: bool = False,
        supervised_mkl: bool = False,
        model_name: str = 'svm',
        mod_options: dict = dict(),
        random_state: int = 0,
        platt_scaling: bool = False,
        bootstrapping: bool = True,
        n_bootstrap_samples: int = 10,
        bootstrap_sample_frac: int = 0.8,
        cv = LeaveOneGroupOut(),
        C: float = 1,
) -> Dict[str, Any]:

    gc_codes = np.asarray(dataset.gc_codes, dtype=object)
    annotations = np.asarray(dataset.annotations, dtype=object)
    groups = dataset.groups
    if hasattr(dataset, 'y'):
        y = dataset.y
    else:
        y = np.asarray([annot.is_cancer for annot in annotations], dtype=int)

    if hasattr(dataset, 'tf'):
        tf = dataset.tf
    else:
        tf = np.copy(dataset.ichor)
        tf[np.logical_and(np.isnan(tf), y == 0)] = np.nanmean(tf[y == 0])

    if y is None:
        if isinstance(annotations[0], Annotation):
            if hasattr(dataset, 'y'):
                y = dataset.y
            else:
                if multi_stage:
                    y = np.asarray([np.nan_to_num(float(annot)) + 1 for annot in annotations], dtype=int)
                    y *= np.asarray([annot.is_cancer for annot in annotations], dtype=int)
                else:
                    y = np.asarray([annot.is_cancer for annot in annotations], dtype=int)
        else:
            if hasattr(dataset, 'y'):
                y = dataset.y
            else:
                y = annotations
            custom_svm = False

    extractor = KernelMatrixExtractor(modality_enums, mod_options=mod_options)
    extractor.preprocess(dataset.data, fit=True)

    n_classes = len(set(y)) + 1

    annotations_target = []
    y_train_pred, y_train_target = [], []
    y_pred, y_target = [], []
    gc_codes_target = []

    #if not bootstrapping:
    #    n_bootstrap_samples = 1
    for iteration in tqdm.tqdm(range(n_bootstrap_samples)):

        if bootstrapping:
            idx = np.random.randint(0, len(dataset), size=int(bootstrap_sample_frac * len(dataset)))
        else:
            idx = np.arange(len(dataset))

        for i, (train_index, test_index) in enumerate(cv.split(np.zeros((len(dataset), 1))[idx], y[idx], groups[idx])):

            try:

                train_index = idx[train_index]
                test_index = idx[test_index]

                y_train = y[train_index]
                y_test = y[test_index]
                tf_train = tf[train_index]

                # Compute one kernel matrix for each modality
                extractor = KernelMatrixExtractor(modality_enums, mod_options=mod_options)
                key_blacklist = {'CpG Methylation (deconvolution, reference)'}
                Ks_train = extractor.fit_transform(
                    Dataset.subset_data(dataset.data, train_index),
                    y=y_train
                )
                Ks_test = extractor.transform(
                    Dataset.subset_data(dataset.data, test_index),
                )

                # Combine the modalities into a single kernel matrix using MKL
                mkl = MultipleKernelLearning(supervised=supervised_mkl)
                K_train = mkl.fit_transform(Ks_train, y_train)
                K_test = mkl.transform(Ks_test)
                assert not np.any(np.isnan(K_train))
                assert not np.any(np.isinf(K_train))

                if model_name == 'stage-svm':
                    model = MultiStageSVM(C=C, C_aux=C, random_state=random_state, platt_scaling=platt_scaling)
                    model.fit(K_train, annotations[train_index])
                elif model_name == 'svm-reimpl':
                    model = SVM(C=C, random_state=random_state, platt_scaling=platt_scaling)
                    model.fit(K_train, y_train)
                elif model_name == 'svm':
                    model = SVC(kernel='precomputed', probability=platt_scaling, C=C, class_weight='balanced', tol=1e-6, random_state=random_state)
                    #model = SVC(kernel='precomputed', class_weight='balanced', probability=True, C=C, tol=1e-6, random_state=random_state)
                    model.fit(K_train, y_train)
                elif model_name == 'carsvm':
                    model = CARSVM(C_c=C, C_r=C*0.05, epsilon=0.1, platt_scaling=platt_scaling)
                    model.fit(K_train, np.asarray([y_train, tf_train]).T)
                else:
                    raise NotImplementedError(f'Unknown model "{model_name}"')

                if platt_scaling:
                    y_train_pred.append(pad_predicted_proba(model.predict_proba(K_train), y_train, n_classes))
                    y_pred.append(pad_predicted_proba(model.predict_proba(K_test), y_train, n_classes))
                else:
                    y_pred_ = model.decision_function(K_train)
                    y_train_pred.append(
                        pad_predicted_proba(
                            np.asarray([-y_pred_, y_pred_]).T,
                            y_train,
                            n_classes
                        )
                    )
                    y_pred_ = model.decision_function(K_test)
                    y_pred.append(
                        pad_predicted_proba(
                            np.asarray([-y_pred_, y_pred_]).T,
                            y_train,
                            n_classes
                        )
                    )

                y_target.append(y_test)
                y_train_target.append(y_train)
                annotations_target.append(annotations[test_index])
                gc_codes_target.append(gc_codes[test_index])
            except Exception as e:
                print(traceback.format_exc())
                continue

    if len(gc_codes_target) > 0:
        gc_codes_target = np.squeeze(np.concatenate(gc_codes_target, axis=0))
        y_train_target = np.squeeze(np.concatenate(y_train_target, axis=0))
        y_train_pred = np.squeeze(np.concatenate(y_train_pred, axis=0))
        y_target = np.squeeze(np.concatenate(y_target, axis=0))
        y_pred = np.squeeze(np.concatenate(y_pred, axis=0))
        annotations_target = np.squeeze(np.concatenate(annotations_target, axis=0))

        if (not multi_stage) and (len(y_train_pred.shape) == 2):
            y_train_pred = np.sum(y_train_pred[:, 1:], axis=1)
        if (not multi_stage) and (len(y_pred.shape) == 2):
            y_pred = np.sum(y_pred[:, 1:], axis=1)

    results = {
        'gc-codes': list(gc_codes_target),
        'annotations_target': annotations_target,
        'y_target': y_target,
        'y_pred': y_pred
    }

    results['roc'] = roc_curve((y_target > 0).astype(int), y_pred)
    results['pr'] = precision_recall_curve((y_target > 0).astype(int), y_pred)
    if len(y_train_pred.shape) == 2:
        results['train-mcc'] = matthews_corrcoef(y_train_target, np.argmax(y_train_pred, axis=1))
    else:
        results['train-mcc'] = matthews_corrcoef(y_train_target, (y_train_pred > 0.5).astype(int))
    if len(y_pred.shape) == 2:
        results['mcc'] = matthews_corrcoef(y_target, np.argmax(y_pred, axis=1))
    else:
        results['mcc'] = matthews_corrcoef(y_target, (y_pred > 0.5).astype(int))

    if multi_stage:
        results['train-auc'] = roc_auc_score(y_train_target, y_train_pred, multi_class='ovr')
        results['auc'] = roc_auc_score(y_target, y_pred, multi_class='ovr')
        results['train-aupr'] = average_precision_score(y_train_target, y_train_pred)
        results['aupr'] = average_precision_score(y_target, y_pred)
    else:
        y_train_target = (y_train_target > 0).astype(int)
        y_target = (y_target > 0).astype(int)
        results['train-auc'] = roc_auc_score(y_train_target, y_train_pred)
        results['auc'] = roc_auc_score(y_target, y_pred)
        results['train-aupr'] = average_precision_score(y_train_target, y_train_pred)
        results['aupr'] = average_precision_score(y_target, y_pred)

        print(f'Train MCC: {results["train-mcc"]}')
        print(f'Train AUROC: {results["train-auc"]}')
        print(f'Train AUPR: {results["train-aupr"]}')
        print(f'MCC: {results["mcc"]}')
        print(f'AUROC: {results["auc"]}')
        print(f'AUPR: {results["aupr"]}')
        print('')

    return results


def mkl_train(
        dataset: Dataset,
        modality_enums: List[Modalities],
        multi_stage: bool = False,
        model_name: str = 'svm-reimpl',
        C: float = 1,
        platt_scaling: bool = False,
        supervised_mkl: bool = False,
        mod_options: dict = dict(),
        random_state: int = 0xCAFE
) -> Tuple[MultipleKernelLearning, Union[SVC, SVM]]:

    annotations = dataset.annotations
    groups = dataset.groups
    y = dataset.y
    tf = dataset.tf

    extractor = KernelMatrixExtractor(modality_enums, mod_options=mod_options)
    Ks = extractor.fit_transform(dataset.data, y)

    mkl = MultipleKernelLearning(supervised=supervised_mkl)
    K = mkl.fit_transform(Ks, y=y)

    weights = np.full(len(K), 1. / len(K))
    weights /= np.sum(weights)
    weights = weights[:, np.newaxis, np.newaxis]
    K = np.sum(weights * K, axis=0)

    if model_name == 'stage-svm':
        model = MultiStageSVM(C=C, C_aux=0.1*C, random_state=random_state, platt_scaling=platt_scaling)
        model.fit(K, annotations)
    elif model_name == 'svm-reimpl':
        model = SVM(C=C, random_state=random_state, platt_scaling=platt_scaling)
        model.fit(K, y)
    elif model_name == 'svm':
        model = SVC(kernel='precomputed', probability=platt_scaling, C=C, class_weight='balanced', tol=1e-6, random_state=random_state)
        #model = SVC(kernel='precomputed', class_weight='balanced', probability=True, C=C, tol=1e-6, random_state=random_state)
        model.fit(K, y)
    elif model_name == 'carsvm':
        model = CARSVM(C_c=C, C_r=C*0.1, epsilon=0.05, platt_scaling=platt_scaling)  # 1, 0.05
        model.fit(K, np.asarray([y, tf]).T)
    else:
        raise NotImplementedError(f'Unknown model "{model_name}"')

    return extractor, mkl, model


def mkl_validate(
        dataset: Dataset,
        mkl: MultipleKernelLearning,
        model: Union[SVM, SVC],
        regression: bool = False,
        multi_stage: bool = False
) -> Dict[str, Any]:

    annotations = dataset.annotations
    groups = dataset.groups

    K = mkl.transform(dataset.data)

    if isinstance(annotations[0], Annotation):
        if hasattr(dataset, 'y'):
            y = dataset.y
        else:
            if multi_stage:
                y = np.asarray([np.nan_to_num(float(annot)) + 1 for annot in annotations], dtype=int)
                y *= np.asarray([annot.is_cancer for annot in annotations], dtype=int)
            else:
                y = np.asarray([annot.is_cancer for annot in annotations], dtype=int)
    else:
        if hasattr(dataset, 'y'):
            y = dataset.y
        else:
            y = annotations

    weights = np.full(len(K), 1. / len(K))
    weights /= np.sum(weights)
    weights = weights[:, np.newaxis, np.newaxis]
    K = np.sum(weights * K, axis=0)

    if isinstance(model, SVM):
        y_pred = model.predict_proba(K)
    else:
        if regression:
            y_pred = model.predict(K)
        else:
            y_pred = model.predict_proba(K)

    y_target = y
    annotations_target = annotations

    results = {
        'annotations_target': annotations_target,
        'y_target': y_target,
        'y_pred': y_pred
    }

    if regression:
        results['r2'] = r2_score(y_target, y_pred)
        print(f'r2: {results["r2"]}')
    else:
        results['roc'] = roc_curve((y_target > 0).astype(int), y_pred[:, 1])
        results['pr'] = precision_recall_curve((y_target > 0).astype(int), y_pred[:, 1])
        results['mcc'] = matthews_corrcoef(y_target, np.argmax(y_pred, axis=1))

        if multi_stage:
            results['auc'] = roc_auc_score(y_target, y_pred, multi_class='ovr')
            results['aupr'] = average_precision_score(y_target, y_pred)
        else:
            results['auc'] = roc_auc_score(y_target, y_pred[:, 1])
            results['aupr'] = average_precision_score(y_target, y_pred[:, 1])

        print(f'MCC: {results["mcc"]}')
        print(f'AUROC: {results["auc"]}')
        print(f'AUPR: {results["aupr"]}')
        print('')

        auc = results['auc']
        roc = results['roc']
        aupr = results['aupr']
        pr = results['pr']

    return results
