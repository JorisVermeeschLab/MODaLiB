from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

from modalib.core import Sample


def parse_mbias(path: str) -> tuple[np.ndarray, ...]:
    """Parse MethylDackel mbias output into OT-R1, OT-R2, OB-R1, OB-R2."""
    df = pd.read_csv(path, sep=r"\s+")
    df["ratio"] = df["nMethylated"] / (df["nMethylated"] + df["nUnmethylated"])
    arrays = []
    for strand, read in [("OT", 1), ("OT", 2), ("OB", 1), ("OB", 2)]:
        subset = df[(df["Strand"] == strand) & (df["Read"] == read)].sort_values("Position")
        arrays.append(subset["ratio"].to_numpy(dtype=float))
    return tuple(arrays)


def infer_mbias_bounds(
    ratios,
    strand="OT",
    central_fraction=0.45,
    window=5,
    min_deviation_5p=0.064,
    min_deviation_3p=0.069,
    mad_factor=1.0,
    min_retained_fraction=0.3,
):
    """Return 0-based, end-exclusive M-bias inclusion bounds."""

    y = np.asarray(ratios, dtype=float)
    n = len(y)

    if n < 2 * window + 2:
        raise ValueError("Read profile too short")

    valid = np.isfinite(y) & (y >= 0) & (y <= 1)
    y = np.where(valid, y, np.nan)

    # Estimate central plateau
    margin = int(n * (1 - central_fraction) / 2)
    central = y[margin:n - margin]

    if np.isfinite(central).mean() < 0.8:
        raise ValueError("Insufficient central coverage")

    baseline = np.nanmedian(central)
    mad = np.nanmedian(np.abs(central - baseline))
    noise = mad_factor * 1.4826 * mad

    threshold_5p = max(min_deviation_5p, noise)
    threshold_3p = max(min_deviation_3p, noise)

    # Rolling median
    smooth = (
        pd.Series(y)
        .rolling(window, center=True, min_periods=1)
        .median()
        .to_numpy()
    )

    # Directional detection
    if strand == "OT":
        bad_5p = (~valid) | (smooth > baseline + threshold_5p)
        bad_3p = (~valid) | (smooth < baseline - threshold_3p)

    elif strand == "OB":
        bad_5p = (~valid) | (smooth < baseline - threshold_5p)
        bad_3p = (~valid) | (smooth > baseline + threshold_3p)

    else:
        raise ValueError(f"Unknown strand: {strand}")

    center = n // 2

    # Identify nearest biased positions
    left_bad = np.flatnonzero(bad_5p[:center])
    right_bad = np.flatnonzero(bad_3p[center:]) + center

    left = left_bad[-1] + 1 if len(left_bad) else 0
    right = right_bad[0] if len(right_bad) else n

    if right <= left:
        raise ValueError("No stable region found")

    if (right - left) / n < min_retained_fraction:
        raise ValueError("Too much of the read would be trimmed")

    if right == len(ratios):
        right = 0

    return int(left), int(right)


def infer_bounds(r1_ot, r2_ot, r1_ob, r2_ob):
    r1_ot_start, r1_ot_end = infer_mbias_bounds(
        r1_ot, strand="OT", min_deviation_5p=0.065, min_deviation_3p=0.069,
        window=5, central_fraction=0.7
    )
    r2_ot_start, r2_ot_end = infer_mbias_bounds(
        r2_ot, strand="OT", min_deviation_5p=0.062, min_deviation_3p=0.069,
        window=5, central_fraction=0.4
    )
    r1_ob_start, r1_ob_end = infer_mbias_bounds(
        r1_ob, strand="OB", min_deviation_5p=0.064, min_deviation_3p=0.066,
        window=5, central_fraction=0.5
    )
    r2_ob_start, r2_ob_end = infer_mbias_bounds(
        r2_ob, strand="OB", min_deviation_5p=0.061, min_deviation_3p=0.071,
        window=3, central_fraction=0.4
    )
    return (
        r1_ot_start, r1_ot_end, r2_ot_start, r2_ot_end,
        r1_ob_start, r1_ob_end, r2_ob_start, r2_ob_end,
    )


# Get methylation inclusion bounds
ot_r1, ot_r2, ob_r1, ob_r2 = parse_mbias(snakemake.input.mbias)
bounds = infer_bounds(ot_r1, ot_r2, ob_r1, ob_r2)

out = Path(snakemake.output.modalib)
out.parent.mkdir(parents=True, exist_ok=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

cmd = [
    str(snakemake.params.tool),
    snakemake.input.reference,
    snakemake.input.bam,
    "-o", str(out),
    "--cconv",
    "--OT", ",".join(map(str, bounds[:4])),
    "--OB", ",".join(map(str, bounds[4:])),
    "-m", str(out) + ".ms.txt",
    "-x", "3"
]

with log_path.open("w") as log:
    log.write("RUNNING: " + repr(cmd) + "\n")
    log.flush()
    subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
if not out.is_file() or out.stat().st_size == 0:
    raise RuntimeError(f"BamConvertor did not create a nonempty .modalib file: {out}")

# Check that file is not corrupted
sample = Sample(str(out))
if sample[len(sample) - 1]["start_chr_id"] != 0:
    raise RuntimeError(f"File {modalib_filepath} is corrupted.")
