#!/bin/bash -l
#
# NAMD on Piz Daint
#
# 1 MPI task per node, 24 OpenMP threads per task with hyperthreading (--ntasks-per-core=2)
#
#SBATCH --job-name="namd"
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
# we do not use the GPU
##SBATCH --constraint=gpu
#SBATCH --constraint=cpu
#========================================
# load modules and run simulation

# cpu modules
module load daint-mc
module load NAMD/2.13-CrayIntel-19.10
# gpu modules
#module load daint-gpu

# instructions
# https://user.cscs.ch/computing/applications/namd/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# energy minimisation
srun namd2 +idlepoll +ppn $[SLURM_CPUS_PER_TASK-1] min.namd > min.log
# equilibriate with different constraints
srun namd2 +idlepoll +ppn $[SLURM_CPUS_PER_TASK-1] eq_step1.namd > eq_step1.log
srun namd2 +idlepoll +ppn $[SLURM_CPUS_PER_TASK-1] eq_step2.namd > eq_step2.log
srun namd2 +idlepoll +ppn $[SLURM_CPUS_PER_TASK-1] eq_step3.namd > eq_step3.log
srun namd2 +idlepoll +ppn $[SLURM_CPUS_PER_TASK-1] eq_step4.namd > eq_step4.log
# run the production
srun namd2 +idlepoll +ppn $[SLURM_CPUS_PER_TASK-1] prod.namd > prod.log