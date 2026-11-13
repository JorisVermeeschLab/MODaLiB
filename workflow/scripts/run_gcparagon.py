from pathlib import Path
import shutil
import subprocess


wc = snakemake.wildcards.sample
out_bam = Path(snakemake.output.bam).resolve()
out_bai = Path(snakemake.output.bai).resolve()
native_out = out_bam.parent / f"{wc}.gcparagon_native"
native_out.mkdir(parents=True, exist_ok=True)
out_bam.parent.mkdir(parents=True, exist_ok=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

cmd = [
    "python",
    str(snakemake.params.tool),
    "--bam", str(Path(snakemake.input.bam).resolve()),
    "--out-dir", str(native_out.resolve()),
    "--output-bam",
    "--two-bit-reference-genome", str(snakemake.input.reference)
]

with log_path.open("w") as log:
    log.write("RUNNING: " + repr(cmd) + "\n")
    log.flush()
    subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)

    # Find the GCparagon output BAM
    candidates = [
        p for p in native_out.rglob("*.GCtagged.bam")
        if p.is_file()
    ]

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one GCtagged BAM under {native_out}, "
            f"found {len(candidates)}: {candidates}"
        )

    native_bam = candidates[0]

    # Find the corresponding BAM index
    native_bai = Path(str(native_bam) + ".bai")

    if not native_bai.exists():
        raise FileNotFoundError(
            f"GCparagon BAM index not found: {native_bai}"
        )

    # Copy BAM to Snakemake output
    shutil.copy2(native_bam, out_bam)

    # Copy BAI to Snakemake output
    shutil.copy2(native_bai, out_bai)

    log.write(f"Copied BAM: {native_bam} -> {out_bam}\n")
    log.write(f"Copied BAI: {native_bai} -> {out_bai}\n")

    # Delete original folder (since we copied the bam and bai files)
    shutil.rmtree(native_out)
