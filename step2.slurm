#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --time=7-00:00:00
#SBATCH --partition=batch_long
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 15
#SBATCH --mem=128G
#SBATCH -A lp_joris_vermeesch
#SBATCH -J mkl
#SBATCH -o mkl.out
#SBATCH -e mkl.err
#SBATCH --mail-user=antoine.passemiers@kuleuven.be
#SBATCH --mail-type=END,FAIL

# Interactive: srun -A lp_joris_vermeesch --cluster=wice --mem=32G --cpus-per-task 16 --time 1:00:00 --pty bash

module --force purge
module use /apps/leuven/skylake/2018a/modules/all
module load Python/3.6.4-foss-2018a
module load matplotlib/2.1.2-foss-2018a-Python-3.6.4
module load GDB/8.1-foss-2018a

# Compilation
cd /lustre1/project/stg_00019/research/Antoine/MODaLiB/step2-fragment-level-analysis
python setup.py install --user

cd ..
python step2.py
