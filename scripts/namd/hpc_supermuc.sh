#!/bin/bash
# 6,336 Thin compute nodes each with 48 cores and 96 GB memory
# 144 Fat compute nodes each 48 cores and 768 GB memory per node
# https://doku.lrz.de/display/PUBLIC/SuperMUC-NG
# docs.google.com/spreadsheets/d/1jbWBUUz3w8_mbxyxJWGOP0znV1gjBBY5DqMHF_RpY7Y/edit#gid=1795899925

#SBATCH --job-name="namd"
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#SBATCH --time=23:20:00
#SBATCH --no-requeue

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48

#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn98ve
#SBATCH --partition=general # test, micro, general, large or fat
#========================================
module load slurm_setup
module load namd

# energy minimisation
mpiexec -n $SLURM_NTASKS namd2 min.namd > min.log
# equilibriate with different constraints
mpiexec -n $SLURM_NTASKS namd2 eq_step1.namd > eq_step1.log
mpiexec -n $SLURM_NTASKS namd2 eq_step2.namd > eq_step2.log
mpiexec -n $SLURM_NTASKS namd2 eq_step3.namd > eq_step3.log
mpiexec -n $SLURM_NTASKS namd2 eq_step4.namd > eq_step4.log
# run the production
mpiexec -n $SLURM_NTASKS namd2 prod.namd > prod.log