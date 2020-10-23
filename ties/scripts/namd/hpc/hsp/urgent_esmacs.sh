#!/bin/bash
# define an array job
#BSUB -J namd[1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 96 # each job run 24 replicas
#BSUB -q compbiomed
##BSUB -q scafellpikeSKL # can be used a litte, 10x4, try 65x1 (does not work). Seems to allow only for 10 at a time.
##BSUB -q scafellpikeKNL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
##BSUB -q universeScafellPike - cannot be used
##BSUB -q scafellpikeI - cannot be used
#BSUB -W 30:00
#BSUB -x

#Load modules
source /etc/profile.d/modules.sh
module load namd-gcc/2.12

# cd $HCBASE/urg_fgesmacs/PLPro-select/6w9c

# from 0 to 2 inclusive
for step in {0..2}; do
    # use a specific drugs
    drug="l2105_rep3_frame_12"
    # each job uses 24 replicas,
    echo "are we in the right place?" `ls`
    # charmrun +p<procs> ++mpiexec ++remote-shell mympiexec namd2 <configfile>

    # PE are worker threads
    # ppn number of PEs per process, can be used only in SMP mode?
    mpiexec -np 3 namd2 +setcpuaffinity +pemap 1-31 +commap 0 +replicas 3 +stdout $drug/rep-%d.log $drug/replica-confs/eq$step-replicas.conf
done
#
#for step in {1..1}; do
#    i=0
#    for drug in `ls -d l*-*`; do
#        ibrun -n 24 -o $((i*24)) namd2 +setcpuaffinity ++ppn 55 +pemap 1-55 +commap 0 +replicas 24 $drug/replica-confs/sim$step-replicas.conf &
#        sleep 2
#        ((i++))
#    done
#    wait
#done