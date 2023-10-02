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

cd /lustre1/project/stg_00019/research/Antoine/MODaLiB

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:step1-bam-processing/dependencies/htslib-1.18/

python step1.py -lb 300
