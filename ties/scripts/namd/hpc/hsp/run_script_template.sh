#!/bin/bash
# define an array job
#BSUB -J anatyk[1-1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 32
#BSUB -q compbiomed
#BSUB -W 15:00
#BSUB -x

mkdir analysis

#Load modules
source /etc/profile.d/modules.sh
module load python3/3.6.2



sample=250

# how many sample replicas
for i in {{1..20}}
do
    python3 ddg_replication.py $sample {transformation} &
done

wait