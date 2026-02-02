#!/bin/bash
#
#SBATCH --job-name=submission
#SBATCH --output=result.txt
#
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --time=120:00

module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

srun -n 48 ./build/numsim_parallel Task4/scenarios/hugeWing4.txt