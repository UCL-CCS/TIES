#!/bin/bash
# define an array job
#BSUB -J namd[1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 32
#BSUB -q compbiomed
##BSUB -q scafellpikeSKL # can be used a litte, 10x4, try 65x1 (does not work). Seems to allow only for 10 at a time.
##BSUB -q scafellpikeKNL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
##BSUB -q universeScafellPike - cannot be used
##BSUB -q scafellpikeI - cannot be used
#BSUB -W 7:00
#BSUB -x



module load python3/3.6.2

for dir in *; do
   if [ -d $dir  ]; then
      cd $dir
          python ../ddg.py > ddg.out &
      cd ..
   fi
done