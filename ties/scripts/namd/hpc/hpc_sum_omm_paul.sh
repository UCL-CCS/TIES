#!/bin/bash
#BSUB -P CHM155_001
#BSUB -W 0:20
#BSUB -nnodes 1
#BSUB -alloc_flags "smt2"
#BSUB -J Test_openmm
#BSUB -o Test_openmm.o%J
#BSUB -e Test_openmm.e%J

export OMP_NUM_THREADS=8
WDIR="$MEMBERWORK/chm155/p/openmm/charmm-gui-9085832103/openmm"
NNODES="1"
NRS="1"
NTASKS="1"
NCORES="7" # Physical Cores per Resource Set
NGPUS="1"

# source ~/.bash_profile
# source ~/.bashrc
# module load cuda/10.1.243
module load gcc/6.4.0  spectrum-mpi/10.3.1.2-20200121
# module load gromacs/2020

nvidia-smi > log_nvidiasmi
env > log_env

echo -e "Nodes\t$NNODES\nMPI tasks\t$NTASKS\nCores\t$NCORES\nGPUs\t$NGPUS\n"
# cd $MEMBERWORK/chm155/miniconda/envs/openmm/share/openmm/examples
# jsrun -n $NRS -a $NTASKS -c $NCORES -g $NGPUS gmx_mpi grompp -f ../../production.mdp -c ../post-anneal.gro -r ../post-anneal.gro -n ../../index.ndx -p ../../topol.top -o production.tpr

jsrun -n $NRS -a $NTASKS -c $NCORES -g $NGPUS gmx_mpi mdrun -v -deffnm production -nb gpu -pme gpu -pmefft gpu -bonded gpu
