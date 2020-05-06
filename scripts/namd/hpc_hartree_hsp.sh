#!/bin/bash
# namd examples: https://stfc.service-now.com/hartreecentre?id=kb_article&sys_id=1c8b4bfcdb440c5034836c16ca961930
# https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_command_ref/bsub.heading_options.1.html

#BSUB -J respties
#BSUB -o std.%J.o
#BSUB -e std.%J.e
#BSUB -R "span[ptile=64]"
# 64 cpus per node? ie two cpus 32 for each? (threads I assume)
# 64 cpu * 65 nodes
# this should be ptile * number of nodes,
#BSUB -n 128
##BSUB -q scafellpikeSKL
#BSUB -q compbiomed
#BSUB -W 1:00
#BSUB -x
##BSUB -cwd # this does not work

# go to the right directionary
cd $HCBASE/resp/tyk2_l5_l16

source /etc/profile.d/modules.sh
module load namd-gcc/2.12
#module load namd/2.13

#	mpiexec.hydra -np 128 namd2 +pemap 1-31 +commap 0 ../eq.conf > eq.log &

# We use timeout to ensure that each NAMD Process is terminated, even if it gets stuck
# 60 seconds * 60 minutes * 4 hours = 4 hours timeout for the ligand

echo "lsb hosts:"
echo $LSB_HOSTS
# count how many processors are allocated
NP=0
for TOKEN in $LSB_HOSTS
do
   ((NP++))
done

# echo $env

TASKS_PER_JOB=32
# use one node for each simulation
NP=64
NODES_PER_JOB=1
BASE_DIR=`pwd`

function schedule_system() {
    local lambda=$1
	local replica=$2

	echo "Scheduling ligand and complex for lambda $lambda replica $replica"
	LIG_PATH="lig/lambda_$lambda/rep$replica"
	COMPLEX_PATH="complex/lambda_$lambda/rep$replica"
	# We use timeout to ensure that each NAMD ligand process is terminated in case it gets stuck
    # 60 seconds * 60 minutes * 4 hours = 4 hours timeout for the ligand
    local LIG_TIMEOUT=$(( 60 * 60 * 4 ))

    # move to the correct ligand directory
    cd $BASE_DIR
    cd $LIG_PATH
	(

	    # first is the ligand part with the timeout
	    ( timeout $LIG_TIMEOUT mpiexec -n $NP namd2    min.namd > min.log &&
        mpiexec -n $NP namd2    eq_step1.namd > eq_step1.log &&
        mpiexec -n $NP namd2    eq_step2.namd > eq_step2.log &&
        mpiexec -n $NP namd2    eq_step3.namd > eq_step3.log &&
        mpiexec -n $NP namd2    eq_step4.namd > eq_step4.log &&
        mpiexec -n $NP namd2    prod.namd > prod.log &&
        echo "Finished ligand: lambda $lambda replica $replica"
        ) ;
        (
        # then the complex part, no timeout
        cd $BASE_DIR && cd $COMPLEX_PATH &&
        mpiexec -n $NP namd2    min.namd > min.log &&
        mpiexec -n $NP namd2    eq_step1.namd > eq_step1.log &&
        mpiexec -n $NP namd2    eq_step2.namd > eq_step2.log &&
        mpiexec -n $NP namd2    eq_step3.namd > eq_step3.log &&
        mpiexec -n $NP namd2    eq_step4.namd > eq_step4.log &&
        mpiexec -n $NP namd2    prod.namd > prod.log &&
        echo "Finished complex: lambda $lambda replica $replica"
        )
    ) &
}

#declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
#declare -a replicas=(1 2 3 4 5)

declare -a lambdas=(0.50)
declare -a replicas=(1 2)

for lambda in "${lambdas[@]}"; do
    for replica in "${replicas[@]}"; do
       schedule_system $lambda $replica
    done
done

wait