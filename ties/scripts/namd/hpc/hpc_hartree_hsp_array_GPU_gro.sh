#!/bin/bash
# define an array job
#BSUB -J testg[1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=16]"
#BSUB -n 16
#BSUB -q scafellpikeGPU # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
#BUSB -gpu "num=1:mode=exclusive_process"
#BSUB -W 00:20
#BSUB -x

SYSTEM=$HCBASE/test/p/

env > last_env
nvidia-smi > smilog

# Load modules
source /etc/profile.d/modules.sh
module load gromacs/2020.1
#source mpivars.sh
#module load intel/latest

gmx_mpi convert-tpr -s step9_1.tpr -extend 1000000 -o step10
mpis=$(( 1 * 8 ))
#mpiexec -np $mpis
gmx_mpi mdrun -v -deffnm step10 # -gpu_id 0 # -ntomp 4 # -rdd 4 -dds 1

