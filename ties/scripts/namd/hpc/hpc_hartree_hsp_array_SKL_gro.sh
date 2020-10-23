#!/bin/bash
# define an array job
#BSUB -J testp[1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 64
#BSUB -q compbiomed
##BSUB -q scafellpikeSKL
#BSUB -W 35:20
#BSUB -x

SYSTEM=$HCBASE/test/p/chol_psm/R1/equilibrate

env > last_env

# Load modules
source /etc/profile.d/modules.sh
module load gromacs/2020.1
source mpivars.sh
#module load intel/latest


csh run_eq
#gmx_mpi convert-tpr -s step9_1.tpr -extend 1000000 -o step10
#mpis=$(( 1 * 8 ))
#
#gx_mpi mdrun -v -deffnm step10 -ntomp 4 # -rdd 4 -dds 1

