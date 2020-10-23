#!/bin/bash
#BSUB -J mmpbsa
#BSUB -o std.%J.o
#BSUB -e std.%J.e
#BSUB -u "bieniekmat@gmail.com"
#BSUB -B # email when executed/dispatched
#BUSB -N # email when the job terminates
#BSUB -R "span[ptile=1]" ###### 32 cores per node
# this should be ptile * number of nodes,
#BSUB -n 32       ##### for 64 cores, total nodes requested 2, 8320 is 128*65
##BSUB -q compbiomed
#BSUB -q scafellpikeSKL
#BSUB -W 0:10
#BSUB -x

ROOT_WORK=/lustre/scafellpike/local/HT03119/mjm06/mxb57-mjm06/software/amber18/bin
cd $ROOT_WORK
NP=32 # cores per simulation

#Load modules
source /etc/profile.d/modules.sh
#module load namd-gcc/2.12
#module load namd/2.13
module load intel/19.3.199
module load intel_mpi/19.3.199
module load python3/3.6.2

python MMPBSA.py.MPI