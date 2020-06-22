#!/bin/bash
# 6,336 Thin compute nodes each with 48 cores and 96 GB memory
# 144 Fat compute nodes each 48 cores and 768 GB memory per node
# https://doku.lrz.de/display/PUBLIC/SuperMUC-NG
# https://doku.lrz.de/display/PUBLIC/Job+Processing+with+SLURM+on+SuperMUC-NG
# https://doku.lrz.de/display/PUBLIC/NAMD
# https://doku.lrz.de/display/PUBLIC/Job+farming+with+SLURM

#SBATCH --job-name="TIESnamd"
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bieniekmat@gmail.com
#SBATCH --no-requeue

#SBATCH --time=20:30:00
#SBATCH --nodes=195
#SBATCH --ntasks-per-node=48

#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn98ve
#SBATCH --partition=general # test, micro, general, large or fat

#constraints are optional
#--constraint="scratch&work"
#========================================
module load slurm_setup
module load namd

TASKS_PER_JOB=48
NODES_PER_JOB=1
BASE_DIR=`pwd`

function schedule_system() {
    local lambda=$1
	local replica=$2

	# note that this is not required, without it will simply be Base
	local transformation=$3
	local THIS_BASE="$BASE_DIR/$transformation"
	echo "The transformation of this directory is: $THIS_BASE"

	echo "Scheduling ligand and complex for lambda $lambda replica $replica"
	LIG_PATH="lig/lambda_$lambda/rep$replica"
	COMPLEX_PATH="complex/lambda_$lambda/rep$replica"
	# We use timeout to ensure that each NAMD ligand process is terminated in case it gets stuck
    # e.g. 60 seconds * 60 minutes * 4 hours = 4 hours timeout for the ligand
    local LIG_TIMEOUT=$(( 60 * 60 * 6 ))

    # move to the correct ligand directory
    cd $THIS_BASE
    cd $LIG_PATH
    cmd_srun_namd="srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2"
	(

	    # first is the ligand part with the timeout
	    ( timeout $LIG_TIMEOUT ${cmd_srun_namd} min.namd > min.log &&
        ${cmd_srun_namd} eq_step1.namd > eq_step1.log &&
        ${cmd_srun_namd} eq_step2.namd > eq_step2.log &&
        ${cmd_srun_namd} eq_step3.namd > eq_step3.log &&
        ${cmd_srun_namd} eq_step4.namd > eq_step4.log &&
        ${cmd_srun_namd} prod.namd > prod.log &&
        echo "Finished ligand: lambda $lambda replica $replica"
        ) ;
        (
        # then the complex part, no timeout
        cd $THIS_BASE && cd $COMPLEX_PATH &&
        ${cmd_srun_namd} min.namd > min.log &&
        ${cmd_srun_namd} eq_step1.namd > eq_step1.log &&
        ${cmd_srun_namd} eq_step2.namd > eq_step2.log &&
        ${cmd_srun_namd} eq_step3.namd > eq_step3.log &&
        ${cmd_srun_namd} eq_step4.namd > eq_step4.log &&
        ${cmd_srun_namd} prod.namd > prod.log &&
        echo "Finished complex: lambda $lambda replica $replica"
        )
    ) &
}

declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
declare -a replicas=(1 2 3 4 5)

declare -a transformations=(l11_l23 l6_l14 l8_l14 l4_l22)
#declare -a transformations=(.)
for base in "${transformations[@]}"; do
    for lambda in "${lambdas[@]}"; do
        for replica in "${replicas[@]}"; do
           schedule_system $lambda $replica $base
        done
    done
done

# wait for all of them to finish
wait