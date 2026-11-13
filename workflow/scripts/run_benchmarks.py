import json
import subprocess
import time
import shutil
import yaml
from datetime import datetime, timezone
from pathlib import Path


sample = str(snakemake.wildcards.sample)
bam = str(Path(snakemake.input.bam).resolve())
reference = str(Path(snakemake.input.reference).resolve())
chrom_sizes = str(Path(snakemake.input.chrom_sizes).resolve())
delfi_bins = str(Path(snakemake.input.delfi_bins).resolve())
wps_sites = str(Path(snakemake.input.wps_sites).resolve())
gap_file = str(Path(snakemake.input.gap_file).resolve())

manifest = Path(snakemake.output.manifest)
output_root = manifest.parent.resolve()
log_path = Path(snakemake.log[0])

output_root.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

threads = str(snakemake.threads)
results = {}


def run_method(method, command, output, log):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    started_at = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()

    log.write(f"[{method}] Running: {' '.join(command)}\n")
    log.flush()

    subprocess.run(
        command,
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(
            f"{method} did not create a nonempty output: {output}"
        )

    results[method] = {
        "started_at": started_at,
        "duration_seconds": time.perf_counter() - start,
        "outputs": [str(output)],
        "command": command,
    }



def run_griffin(log):
    """
    Run official GRIFFIN GC correction and nucleosome
    profiling workflows for one BAM sample.
    """

    griffin_root = Path(
        snakemake.params.griffin_root
    ).resolve()

    griffin_ref = griffin_root / "Ref"
    griffin_scripts = griffin_root / "scripts"

    griffin_output = output_root / "griffin"
    griffin_output.mkdir(parents=True, exist_ok=True)

    started_at = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()

    # ---------------------------------------------
    # GRIFFIN GC correction
    # ---------------------------------------------

    gc_source = (
        griffin_root / "snakemakes" /
        "griffin_GC_and_mappability_correction"
    )

    gc_work = griffin_output / "gc_workflow"
    shutil.copytree(gc_source, gc_work, dirs_exist_ok=True)

    gc_config_path = gc_work / "config" / "config.yaml"
    gc_samples_path = gc_work / "config" / "samples.yaml"
    gc_cluster_path = gc_work / "config" / "cluster_slurm.yaml"

    gc_results = (griffin_output / "gc_results").resolve()

    with gc_config_path.open() as f:
        gc_config = yaml.safe_load(f)

    gc_config.update({
        "griffin_scripts_dir": str(griffin_scripts),
        "reference_genome": reference,
        "chrom_sizes": chrom_sizes,
        "mappability_bw": str(
            griffin_ref / "k100.Umap.MultiTrackMappability.bw"
        ),
        "encode_exclude": str(
            griffin_ref / "encode_unified_GRCh38_exclusion_list.bed"
        ),
        "centromeres": str(
            griffin_ref / "hg38_centromeres.bed"
        ),
        "gaps": str(
            griffin_ref / "hg38_gaps.bed"
        ),
        "patches": str(
            griffin_ref / "hg38_fix_patches.bed"
        ),
        "alternative_haplotypes": str(
            griffin_ref / "hg38_alternative_haplotypes.bed"
        ),
        "mappable_regions": str(
            griffin_ref /
            "k100_minus_exclusion_lists.mappable_regions.hg38.bed"
        ),
        "genome_GC_frequency": str(
            griffin_ref / "genome_GC_frequency"
        ),
        "mappability_correction": False,
        "out_dir": str(gc_results),
    })

    with gc_config_path.open("w") as f:
        yaml.safe_dump(gc_config, f)

    with gc_samples_path.open("w") as f:
        yaml.safe_dump(
            {"samples": {sample: bam}},
            f
        )

    with gc_cluster_path.open() as f:
        gc_cluster = yaml.safe_load(f)

    gc_cluster["GC_counts"]["ncpus"] = int(threads)

    with gc_cluster_path.open("w") as f:
        yaml.safe_dump(gc_cluster, f)

    gc_command = [
        "snakemake",
        "--snakefile",
        "griffin_GC_and_mappability_correction.snakefile",
        "--cores", "1",
        "--rerun-incomplete",
    ]

    log.write(f"[griffin GC] {gc_command!r}\n")
    log.flush()

    subprocess.run(
        gc_command,
        cwd=gc_work,
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

    gc_bias = (
        gc_results / "GC_bias" /
        f"{sample}.GC_bias.txt"
    )

    if not gc_bias.is_file() or gc_bias.stat().st_size == 0:
        raise RuntimeError(
            f"GRIFFIN GC correction failed: {gc_bias}"
        )

    # ---------------------------------------------
    # GRIFFIN nucleosome profiling
    # ---------------------------------------------

    np_source = (
        griffin_root / "snakemakes" /
        "griffin_nucleosome_profiling"
    )

    np_work = griffin_output / "np_workflow"
    shutil.copytree(np_source, np_work, dirs_exist_ok=True)

    np_config_path = np_work / "config" / "config.yaml"
    np_samples_path = np_work / "config" / "samples.GC.yaml"
    np_sites_path = np_work / "config" / "sites.yaml"
    np_cluster_path = np_work / "config" / "cluster_slurm.yaml"

    np_results = (griffin_output / "np_results").resolve()
    np_tmp = (griffin_output / "tmp").resolve()

    with np_config_path.open() as f:
        np_config = yaml.safe_load(f)

    np_config.update({
        "griffin_scripts_dir": str(griffin_scripts),
        "reference_genome": reference,
        "chrom_sizes_path": chrom_sizes,
        "mappability_bw": str(
            griffin_ref / "k100.Umap.MultiTrackMappability.bw"
        ),
        "encode_exclude": str(
            griffin_ref / "encode_unified_GRCh38_exclusion_list.bed"
        ),
        "centromeres": str(
            griffin_ref / "hg38_centromeres.bed"
        ),
        "gaps": str(
            griffin_ref / "hg38_gaps.bed"
        ),
        "patches": str(
            griffin_ref / "hg38_fix_patches.bed"
        ),
        "alternative_haplotypes": str(
            griffin_ref / "hg38_alternative_haplotypes.bed"
        ),
        "results_dir": str(np_results),
        "tmp_dir": str(np_tmp),
        "samples_yaml": str(np_samples_path),
        "sites_yaml": str(np_sites_path),
        "cluster_slurm": str(np_cluster_path),
    })

    with np_config_path.open("w") as f:
        yaml.safe_dump(np_config, f)

    with np_samples_path.open("w") as f:
        yaml.safe_dump({
            "samples": {
                sample: {
                    "bam": bam,
                    "GC_bias": str(gc_bias),
                }
            }
        }, f)

    with np_sites_path.open("w") as f:
        yaml.safe_dump({
            "site_lists": {
                "CTCF": str(
                    Path(snakemake.input.griffin_sites).resolve()
                )
            }
        }, f)

    with np_cluster_path.open() as f:
        np_cluster = yaml.safe_load(f)

    np_cluster["calc_cov"]["ncpus"] = int(threads)
    np_cluster["merge_sites"]["ncpus"] = int(threads)

    with np_cluster_path.open("w") as f:
        yaml.safe_dump(np_cluster, f)

    np_command = [
        "snakemake",
        "--snakefile",
        "griffin_nucleosome_profiling.snakefile",
        "--cores", "1",
        "--rerun-incomplete",
    ]

    log.write(f"[griffin profiling] {np_command!r}\n")
    log.flush()

    subprocess.run(
        np_command,
        cwd=np_work,
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

    corrected = (
        np_results / sample /
        f"{sample}.GC_corrected.coverage.tsv"
    )
    uncorrected = (
        np_results / sample /
        f"{sample}.uncorrected.coverage.tsv"
    )

    for path in (corrected, uncorrected):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing GRIFFIN output: {path}"
            )

    results["griffin"] = {
        "started_at": started_at,
        "duration_seconds": time.perf_counter() - start,
        "outputs": [
            str(gc_bias),
            str(corrected),
            str(uncorrected),
        ],
        "commands": [gc_command, np_command],
    }


with log_path.open("w", encoding="utf-8") as log:

    # DELFI
    delfi_output = output_root / "delfi" / "delfi.tsv"
    run_method(
        "delfi",
        [
            "finaletoolkit", "delfi",
            bam,
            chrom_sizes,
            reference,
            delfi_bins,
            "--gap-file", gap_file,
            "--no-remove-nocov",
            "--threads", threads,
            "--output", str(delfi_output),
            "--verbose",
        ],
        delfi_output,
        log,
    )

    # Window Protection Score
    wps_output = output_root / "wps" / "wps.bw"
    run_method(
        "wps",
        [
            "finaletoolkit", "wps",
            bam,
            wps_sites,
            "--chrom-sizes", chrom_sizes,
            "--reference", reference,
            "--threads", threads,
            "--output", str(wps_output),
            "--verbose",
        ],
        wps_output,
        log,
    )

    # GRIFFIN
    #run_griffin(log)

# Write manifest atomically
manifest.parent.mkdir(parents=True, exist_ok=True)
manifest_tmp = manifest.with_name(manifest.name + ".tmp")

manifest_tmp.write_text(
    json.dumps(
        {"sample": sample, "methods": results},
        indent=2,
    ) + "\n",
    encoding="utf-8",
)

manifest_tmp.replace(manifest)
