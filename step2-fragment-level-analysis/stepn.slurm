#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A lp_joris_vermeesch
#SBATCH -J nuc
#SBATCH -o nuc.out
#SBATCH -e nuc.err
#SBATCH --mail-user=antoine.passemiers@kuleuven.be
#SBATCH --mail-type=END,FAIL

module --force purge
module use /apps/leuven/skylake/2018a/modules/all
module load Python/3.6.4-foss-2018a
module load matplotlib/2.1.2-foss-2018a-Python-3.6.4
module load GDB/8.1-foss-2018a

cd /lustre1/project/stg_00019/research/Antoine/RECEDE/step2-building-features

python setup.py install --user

# python tmp.py
python create-nucleosome-map.py
