#!/bin/bash
# define an array job
#BSUB -J testg[1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=8]"
#BSUB -n 8
#BSUB -q scafellpikeKNL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
#BSUB -W 00:20
#BSUB -x

SYSTEM=$HCBASE/test/p/

env > last_env

# Load modules
source /etc/profile.d/modules.sh
module load gromacs-knl/2020.1
#source mpivars.sh
#module load intel/latest

gmx_mpi convert-tpr -s step9_1.tpr -extend 1000000 -o step10
mpis=$(( 1 * 8 ))
mpiexec -np $mpis gmx_mpi mdrun -v -deffnm step10 -ntomp 4 # -rdd 4 -dds 1


# 20 ns/day with 64 omp and one node,
