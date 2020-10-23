#!/bin/bash
# define an array job
#BSUB -J analysis[1-1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=64]"
#BSUB -n 64
#BSUB -q compbiomed
##BSUB -q scafellpikeSKL # can be used a litte, 10x4, try 65x1 (does not work). Seems to allow only for 10 at a time.
##BSUB -q scafellpikeKNL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
##BSUB -q universeScafellPike - cannot be used
##BSUB -q scafellpikeI - cannot be used
#BSUB -W 30:00
#BSUB -x

SYSTEM=$HCBASE/resp_ff/
NP=32

module load python3/3.6.2

for D in */;
do
    cd $D
    echo "Next Dir: $D"
#    if ! test -f "ddg.out"; then
        python ~/ddg.py > ddg.out &
#    fi
    cd ..
done


wait