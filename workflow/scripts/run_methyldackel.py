import subprocess
from pathlib import Path


bam = snakemake.input.bam
ref = snakemake.input.ref
output = Path(snakemake.output.mbias)
tool = snakemake.params.tool

output.parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)

cmd = [
    tool,
    "mbias",
    "--txt",
    "--noSVG",
    "-@",
    str(snakemake.threads),
    ref,
    bam,
]

with open(output, "w") as out, open(snakemake.log[0], "w") as log:
    subprocess.run(
        cmd,
        stdout=out,
        stderr=log,
        check=True
    )
