from pathlib import Path
import os
import numpy as np
import scanpy


ROOT = Path(__file__).resolve().parent
REGIONS_FOLDER = str((ROOT / "../../regions").resolve())

out = Path(snakemake.output.npz)
out.parent.mkdir(parents=True, exist_ok=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with log_path.open("w") as log:

    ann = scanpy.read_h5ad(str(snakemake.input.h5ad))

    cell_types = ann.obs.cell_type.to_numpy()
    unique_cell_types = np.unique(cell_types)
    ensembl_ids = ann.var.index.to_numpy()

    means, stds = [], []
    for cell_type in unique_cell_types:
        mask = (cell_types == cell_type)
        mean = np.squeeze(np.asarray(ann.X[mask, :].mean(axis=0)))
        means.append(mean)
        std = np.sqrt(np.mean(np.square(ann.X[mask, :].todense() - mean[np.newaxis, :]), axis=0))
        stds.append(std)
    means = np.asarray(means)
    stds = np.asarray(stds)

    mapping = {}
    with open(os.path.join(str(REGIONS_FOLDER), 'ensembl-ids.txt'), 'r') as f:
        for line in f:
            elements = line.rstrip().split()
            if len(elements) == 1:
                mapping[elements[0]] = ''
            elif len(elements) == 2:
                mapping[elements[0]] = elements[1]

    gene_symbols = np.asarray([mapping[ensembl_id] for ensembl_id in ensembl_ids], dtype=object)

    np.savez(
        str(out),
        mu=means,
        sigma=stds,
        ensembl_ids=ensembl_ids,
        cell_types=unique_cell_types,
        gene_symbols=gene_symbols
    )
