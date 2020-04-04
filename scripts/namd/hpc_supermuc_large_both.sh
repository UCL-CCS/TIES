#!/bin/bash
# 6,336 Thin compute nodes each with 48 cores and 96 GB memory
# 144 Fat compute nodes each 48 cores and 768 GB memory per node
# https://doku.lrz.de/display/PUBLIC/SuperMUC-NG
# https://doku.lrz.de/display/PUBLIC/Job+Processing+with+SLURM+on+SuperMUC-NG
# https://doku.lrz.de/display/PUBLIC/NAMD
# https://doku.lrz.de/display/PUBLIC/Job+farming+with+SLURM

#SBATCH --job-name="namd"
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
#SBATCH --nodes=65
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


function schedule_system() {
	lambda=$1
	replica=$2
	echo "Scheduling ligand and complex for lambda $lambda replica $replica"
	LIG_PATH="lig/lambda_%lambda/rep%replica"
	COMPLEX_PATH="complex/lambda_%lambda/rep%replica"
	# We use timeout to ensure that each NAMD ligand process is terminated in case it gets stuck
    # 60 seconds * 60 minutes * 4 hours = 4 hours timeout for the ligand
    local LIG_TIMEOUT=$(( 60 * 60 * 4 ))
	(

	    # first is the ligand part with the timeout
	    (timeout -v $LOG_TIMEOUT srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $LIG_PATH/min.namd > $LIG_PATH/min.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $LIG_PATH/eq_step1.namd > $LIG_PATH/eq_step1.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $LIG_PATH/eq_step2.namd > $LIG_PATH/eq_step2.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $LIG_PATH/eq_step3.namd > $LIG_PATH/eq_step3.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $LIG_PATH/eq_step4.namd > $LIG_PATH/eq_step4.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $LIG_PATH/prod.namd > $LIG_PATH/prod.log &&
        echo "Finished ligand: lambda $lambda replica $replica"
        ) ;
        (
        # then the complex part, no timeout
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $COMPLEX_PATH/min.namd > $COMPLEX_PATH/min.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $COMPLEX_PATH/eq_step1.namd > $COMPLEX_PATH/eq_step1.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $COMPLEX_PATH/eq_step2.namd > $COMPLEX_PATH/eq_step2.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $COMPLEX_PATH/eq_step3.namd > $COMPLEX_PATH/eq_step3.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $COMPLEX_PATH/eq_step4.namd > $COMPLEX_PATH/eq_step4.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 $COMPLEX_PATH/prod.namd > $COMPLEX_PATH/prod.log &&
        echo "Finished complex: lambda $lambda replica $replica"
        )
    ) &
}

declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
declare -a replicas=(1 2 3 4 5)

for lambda in "${lambdas[@]}"; do
    for replica in "${replicas[@]}"; do
       schedule_system $lambda $replica
    done
done

# wait for all of them to finish
wait