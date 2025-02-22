from typing import Dict, List, Tuple

import numpy as np
from sklearn.preprocessing import StandardScaler

from .deconvolution.da import interdomain_deconvolution


def find_common_genes(
        gene_names: List[str],
        ref_gene_names: List[str]
) -> Tuple[np.ndarray, np.ndarray]:
    ref_gene_dict = {gene_name: i for i, gene_name in enumerate(ref_gene_names)}
    shared_gene_names, idx1, idx2 = [], [], []
    for k, gene_name in enumerate(gene_names):
        if isinstance(gene_name, str) and (gene_name in ref_gene_dict):
            idx1.append(k)
            idx2.append(ref_gene_dict[gene_name])
            shared_gene_names.append(gene_name)
    idx1 = np.asarray(idx1, dtype=int)
    idx2 = np.asarray(idx2, dtype=int)
    return idx1, idx2


def infer_expression(
        data: Dict[str, np.ndarray],
        gene_names: List[str],
        ref_exp: np.ndarray,
        ref_gene_names: List[str],
        model = None
) -> Tuple[np.ndarray, List[str]]:

    #ref_exp = np.log(ref_exp + 0.001)
    ref_exp = np.sqrt(ref_exp)

    idx1, idx2 = find_common_genes(gene_names, ref_gene_names)

    shared_gene_names = np.asarray(gene_names, dtype=object)[idx1]

    X = np.asarray([
        #data['Fragment length (first moment)'] / 300.0 - 0.5,
        #data['CpG Methylation (first moment)'] - 0.7,
        data['Proportion of long fragments'] - 0.5,
        data['Nucleosome positioning'] - 0.5,
        #data['GC-content (first moment)'] - 0.5,
    ]).T
    X = np.squeeze(X)

    print("SHAPE", data['Fragment length (first moment)'])

    X = StandardScaler().fit_transform(X)
    
    alpha, model = interdomain_deconvolution(X[idx1, :], ref_exp[:, idx2], model=model)

    return alpha, shared_gene_names, model
