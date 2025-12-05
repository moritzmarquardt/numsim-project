#!/bin/bash
#
#SBATCH --job-name=submission
#SBATCH --output=result.txt
#
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=3
#SBATCH --time=10:00

module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

srun -n 9 ./build/numsim_parallel lid_driven_cavity.txt