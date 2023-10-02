#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH -A lp_joris_vermeesch
#SBATCH -J extr
#SBATCH -o extr.out
#SBATCH -e extr.err
#SBATCH --mail-user=antoine.passemiers@kuleuven.be
#SBATCH --mail-type=END,FAIL

cd /lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing/

module --force purge
module use /apps/leuven/skylake/2018a/modules/all
module load foss/2018a
module load Python/3.6.4-foss-2018a
module load SAMtools/1.9-foss-2018a

cd /lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing/dependencies/htslib-1.18/
make clean
make

cd /lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing
make clean
make

# zcat file.fast.gz | less
# samtools view GC114135.bam | less -S
# sambamba view -f bam -F "not secondary_alignment and proper_pair and mapping_quality >= 50" -t 20 --subsampling-seed=37 -s 0.1 input.bam -o output.bam
# sort -V -k1,1 -k2,2 /lustre1/project/stg_00019/Methylation_atlas_replicates/raw_bed/lung1.2bp.bed > /lustre1/project/stg_00019/research/Antoine/DeconMeth-Extractor/out/lung1.2bp.bed

cd /lustre1/project/stg_00019/research/Antoine/MODaLiB

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:step1-bam-processing/dependencies/htslib-1.18/

python step1.py -lb 0 -ub 100
