from typing import List, Tuple

import numpy as np
import pandas as pd
import psmpy


def match_samples_deprecated(sample_names: List[str], group: np.ndarray, *confounders: np.ndarray) -> List[str]:

    import rpy2
    from rpy2.robjects.packages import importr
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri as rpyn
    import rpy2.robjects.pandas2ri as pandas2ri

    n_pos_before = int(np.sum(group == 1))
    n_neg_before = int(np.sum(group == 0))

    # Convert data to numpy arrays
    confounders = list(confounders)
    sample_names = np.asarray(sample_names, dtype=object)
    group = np.asarray(group, dtype=float)
    for k in range(len(confounders)):
        confounders[k] = np.asarray(confounders[k], dtype=float)
        assert len(group) == len(confounders[k])
        assert len(sample_names) == len(confounders[k])
        invalid = np.isnan(confounders[k])
        if np.any(invalid):
            print(f'Missing values encountered for {np.sum(invalid)} samples out of {len(sample_names)}')
    assert not np.any(np.isnan(group))

    # Create dataset from groups and confounders
    confounder_names = {f'c{i + 1}' for i in range(len(confounders))}
    df = pd.DataFrame({'group': group, **{c_name: c for c_name, c in zip(confounder_names, confounders)}})
    with robjects.conversion.localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)

    # Data imputation
    robjects.packages.importr('mice')
    mice_func = robjects.r['mice']
    mids = mice_func(r_df)

    # Case/control matching
    robjects.packages.importr('MatchThem')
    matchthem_func = robjects.r['matchthem']
    res = matchthem_func(
        formula=robjects.Formula(f'group ~ {" + ".join(confounder_names)}'),
        datasets=mids,
        method='optimal',
        distance='glm',
        discard='control',
        replace=False,
        ratio=1,
        reestimate=True
    )
    res = res[2][0]

    # Get selected subset
    weights = np.asarray(res[2])#[0, :]
    mask = (weights > 0)
    sample_names = sample_names[mask]
    group = group[mask]
    for i in range(len(confounders)):
        confounders[i] = confounders[i][mask]

    print(f'Controls: {n_neg_before} -> {int(np.sum(group == 0))}')
    print(f'Cases   : {n_pos_before} -> {int(np.sum(group == 1))}')

    return list(sample_names)


def match_samples(sample_names: List[str], group: np.ndarray, *confounders: np.ndarray) -> List[str]:

    n_pos_before = int(np.sum(group == 1))
    n_neg_before = int(np.sum(group == 0))

    # Convert data to numpy arrays
    confounders = list(confounders)
    sample_names = np.asarray(sample_names, dtype=object)
    group = np.asarray(group, dtype=float)
    for k in range(len(confounders)):
        confounders[k] = np.asarray(confounders[k], dtype=float)
        assert len(group) == len(confounders[k])
        assert len(sample_names) == len(confounders[k])
        invalid = np.isnan(confounders[k])
        if np.any(invalid):
            print(f'Missing values encountered for {np.sum(invalid)} samples out of {len(sample_names)}')
    assert not np.any(np.isnan(group))

    data = {
        'gc-code': sample_names,
        'group': group
    }
    for k in range(len(confounders)):
        data[f'C{k + 1}'] = confounders[k]
    df = pd.DataFrame(data)
    psm = psmpy.PsmPy(df, treatment='group', indx='gc-code')
    psm.logistic_ps(balance=True)
    psm.knn_matched(matcher='propensity_logit', replacement=False, caliper=None, drop_unmatched=True)
    
    gc_codes = list(psm.matched_ids['matched_ID']) + list(psm.matched_ids['gc-code'])

    print(gc_codes)

    return gc_codes


def match_and_split_train_test(sample_names: List[str], group: np.ndarray, batch: np.ndarray, *confounders: np.ndarray) -> Tuple[List[str], List[str]]:
    raise NotImplementedError()

