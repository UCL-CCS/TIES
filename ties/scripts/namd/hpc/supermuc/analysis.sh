#!/bin/bash
# 6,336 Thin compute nodes each with 48 cores and 96 GB memory
# 144 Fat compute nodes each 48 cores and 768 GB memory per node
# https://doku.lrz.de/display/PUBLIC/SuperMUC-NG
# https://doku.lrz.de/display/PUBLIC/Job+Processing+with+SLURM+on+SuperMUC-NG
# https://doku.lrz.de/display/PUBLIC/NAMD
# https://doku.lrz.de/display/PUBLIC/Job+farming+with+SLURM

#SBATCH --job-name="TIESanalysis"
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bieniekmat@gmail.com
#SBATCH --no-requeue

#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48

#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn98ve
#SBATCH --partition=test # test, micro, general, large or fat

#constraints are optional
#--constraint="scratch&work"
#========================================
module load slurm_setup
module load python/3.6_intel

for D in */;
do
    cd $D
    echo "Next Dir: $D"
#    if ! test -f "ddg.out"; then
        python ../ddg.py > ddg.out &
#    fi
    cd ..
done

wait