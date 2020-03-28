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

#SBATCH --time=10:30:00
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

# We use timeout to ensure that each NAMD Process is terminated, even if it gets stuck
# 60 seconds * 60 minutes * 4 hours = 4 hours timeout for the ligand
LIG_TIMEOUT=$(( 60 * 60 * 4 ))

function schedule_system() {
    echo "system $1"
	echo "Lambda $2"
	echo "Replica $3"
	echo "Scheduling $1 lambda $2 replica $3"
	(
	    timeout -v $LOG_TIMEOUT srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 min.namd > min.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 eq_step1.namd > eq_step1.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 eq_step2.namd > eq_step2.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 eq_step3.namd > eq_step3.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 eq_step4.namd > eq_step4.log &&
        srun -N $NODES_PER_JOB -n $TASKS_PER_JOB namd2 prod.namd > prod.log &&
        echo "Finished protocol for $1 lambda $2 replica $3"
    ) &
}

# schedule the ligands
ROOT_DIR=complex
cd $ROOT_DIR
    for L in lambda*/ ; do
        cd $L
            for R in rep*/ ; do
                cd $R
                    schedule_system $ROOT_DIR $L $R
                cd ..
            done
        cd ..
    done
cd ..

# wait for all of them to finish
wait