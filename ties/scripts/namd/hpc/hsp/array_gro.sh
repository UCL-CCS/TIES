#!/bin/bash
# define an array job
#BSUB -J c10g
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=16]"
#BSUB -n 16
#BSUB -q scafellpikeGPU
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 48:00
#BSUB -x

SYSTEM=$HCBASE/test/m/c10_1.9nmWaterEq

# Load modules
source /etc/profile.d/modules.sh
module load gromacs-gpu/2020.1

#gmx_mpi grompp -f production.mdp -o production.tpr -c ../anneal/anneal.gro -r ../anneal/anneal.gro -n ../../index -p ../../topol.top
#gmx_mpi mdrun -v -deffnm production -pin on -nb gpu -pme gpu -pmefft gpu -bonded gpu #-nt 16 -ntomp 4 -rdd 4 -dds 1
#gmx_mpi convert-tpr -s nvt5ns.tpr -extend 50000 -o nvt50ns
gmx_mpi mdrun -v -deffnm nvt50ns

wait