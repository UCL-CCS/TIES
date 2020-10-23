#!/bin/bash
#BSUB -J hybrid
#BSUB -o std.%J.o
#BSUB -e std.%J.e
#BSUB -R "span[ptile=32]" ###### 32 cores per node
# this should be ptile * number of nodes,
#BSUB -n 32       ##### for 64 cores, total nodes requested 2, 8320 is 128*65
#BSUB -q compbiomed
#BSUB -W 23:00
#BSUB -x

/lustre/scafellpike/local/HT03119/mjm06/mxb57-mjm06/software/NAMD_2.14_Linux-x86_64-multicore/namd2 +p 32  eq_step1.namd > eq_step1.log

