#!/bin/bash
# define an array job
#BSUB -J testg[1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 32
#BSUB -q scafellpikeSKL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
#BSUB -W 00:20
#BSUB -x

SYSTEM=$HCBASE/test/p/post-anneal/R1

env > last_env
nvidia-smi > smilog

# Load modules
source /etc/profile.d/modules.sh
module load gromacs/2020.1
source mpivars.sh
#module load intel/latest

gmx_mpi grompp -f ../../production.mdp -c ../post-anneal.gro -r ../post-anneal.gro -n ../../index.ndx -p ../../topol.top -o production.tpr

mpiexec.hydra -np 32 gmx_mpi mdrun -v -deffnm production -ntomp 1

