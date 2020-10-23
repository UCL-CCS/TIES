#!/bin/bash
#BSUB -J perf
#BSUB -o std.%J.o
#BSUB -e std.%J.e
#BSUB -u "bieniekmat@gmail.com"
#BSUB -B # email when executed/dispatched
#BUSB -N # email when the job terminates
#BSUB -R "span[ptile=1]" ###### 32 cores per node
# this should be ptile * number of nodes,
#BSUB -n 64       ##### for 64 cores, total nodes requested 2, 8320 is 128*65
#BSUB -q compbiomed
#BSUB -W 1:00
#BSUB -x

ROOT_WORK=/lustre/scafellpike/local/HT03119/mjm06/mxb57-mjm06/testsystem/rep1
cd $ROOT_WORK
NP=64 # cores per simulation

#Load modules
source /etc/profile.d/modules.sh
#module load namd-gcc/2.12
module load namd/2.13

cmd_mpinamd="charmrun +n 64 ++mpiexec ++remote-shell mpiexec \
  /lustre/scafellpike/local/apps/intel/namd/2.13/bin/namd2 \
  ++ppn 30 +pemap 1-15,17-31 +commap 0,16"

# run the next simulation
#cmd_mpinamd="mpiexec -n $NP namd2 +pemap 1-31 +commap 0"

${cmd_mpinamd} min.namd > min.log
${cmd_mpinamd} eq_step1.namd > eq_step1.log
${cmd_mpinamd} eq_step2.namd > eq_step2.log
${cmd_mpinamd} eq_step3.namd > eq_step3.log
${cmd_mpinamd} eq_step4.namd > eq_step4.log
${cmd_mpinamd} prod.namd > prod.log
