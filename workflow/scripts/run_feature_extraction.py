from typing import Optional, Tuple
import subprocess
from pathlib import Path
import os
import numpy as np
import pandas as pd

from modalib.core import extract_features


ROOT = Path(__file__).resolve().parent
REGIONS_FOLDER = str((ROOT / "../../regions").resolve())


def parse_mbias(path: str) -> tuple[np.ndarray, ...]:
    """Parse MethylDackel mbias output into OT-R1, OT-R2, OB-R1, OB-R2."""
    df = pd.read_csv(path, sep=r"\s+")
    df["ratio"] = df["nMethylated"] / (df["nMethylated"] + df["nUnmethylated"])
    arrays = []
    for strand, read in [("OT", 1), ("OT", 2), ("OB", 1), ("OB", 2)]:
        subset = df[(df["Strand"] == strand) & (df["Read"] == read)].sort_values("Position")
        arrays.append(subset["ratio"].to_numpy(dtype=float))
    return tuple(arrays)


def process_modalib_file(
    in_filepath: str,
    in_jag_filepath: str,
    out_filepath: str,
    region_list_filepath: str,
    mbias_curves: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    target_n_fragments: Optional[int] = None,
    downsample: bool = False
) -> None:
    kwargs = {}
    if downsample:
        kwargs["target_n_fragments"] = target_n_fragments
    nuc_filepath = str(snakemake.params.nuc_reference)
    matrices = extract_features([in_filepath], nuc_filepath, region_list_filepath, weighting=True, **kwargs)
    for key in matrices:
        matrices[key] = np.squeeze(matrices[key])
    matrices['OT-R1'] = mbias_curves[0]
    matrices['OT-R2'] = mbias_curves[1]
    matrices['OB-R1'] = mbias_curves[2]
    matrices['OB-R2'] = mbias_curves[3]

    # Add methylation info
    M = np.zeros((501, 501), dtype=float)
    U = np.zeros((501, 501), dtype=float)
    with open(in_jag_filepath, 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            elements = line.rstrip().split(',')
            if len(elements) == 4:
                i = int(elements[0])
                j = int(elements[1])
                M[i, j] = int(elements[2])
                U[i, j] = int(elements[3])
    matrices['Overall-methylated-counts'] = M
    matrices['Overall-unmethylated-counts'] = U

    np.savez(out_filepath, **matrices)


modalib_filepath = str(snakemake.input.modalib)
modalib_jag_filepath = str(snakemake.input.modalib_jag)

out = Path(snakemake.output.windows).parent
out.mkdir(parents=True, exist_ok=True)

Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)

mbias_curves = parse_mbias(str(snakemake.input.mbias))

kwargs = {
    "downsample": snakemake.params.downsample,
    "target_n_fragments": snakemake.params.target_n_fragments
}

with open(snakemake.log[0], "w") as log:

    process_modalib_file(
        modalib_filepath,
        modalib_jag_filepath,
        os.path.join(str(out), f'ctcf-noXY.npz'),
        os.path.join(REGIONS_FOLDER, 'ctcf-no-XY.csv'),
        mbias_curves,
        **kwargs
    )

    #process_modalib_file(
    #    modalib_filepath,
    #    modalib_jag_filepath,
    #    os.path.join(str(out), f'50kb.npz'),
    #    os.path.join(REGIONS_FOLDER, 'all-regions.csv'),
    #    mbias_curves,
    #    **kwargs
    #)

    process_modalib_file(
        modalib_filepath,
        modalib_jag_filepath,
        os.path.join(str(out), f'50kb-noXY.npz'),
        os.path.join(REGIONS_FOLDER, 'all-regions-no-XY.csv'),
        mbias_curves,
        **kwargs
    )

    for key in ['promoters', 'tss', 'genebodies']:
        process_modalib_file(
            modalib_filepath,
            modalib_jag_filepath,
            os.path.join(str(out), f'{key}-noXY.npz'),
            os.path.join(REGIONS_FOLDER, f'{key}-no-XY.csv'),
            mbias_curves,
            **kwargs
        )

    process_modalib_file(
        modalib_filepath,
        modalib_jag_filepath,
        os.path.join(str(out), f'metdecode-noXY.npz'),
        os.path.join(REGIONS_FOLDER, 'metdecode', 'dmrs.csv'),
        mbias_curves,
        **kwargs
    )

    process_modalib_file(
        modalib_filepath,
        modalib_jag_filepath,
        os.path.join(str(out), f'meth.npz'),
        os.path.join(REGIONS_FOLDER, 'metdecode', 'dmrs-methmarkerdb-colorectal.csv'),
        mbias_curves,
        **kwargs
    )
