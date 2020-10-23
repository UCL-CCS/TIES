#!/bin/bash
# define an array job
#BSUB -J namd[1-65]
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

NP=32

#Load modules
source /etc/profile.d/modules.sh
module load namd-gcc/2.12

# look up the index for this job instance
if [ -z $LSB_REMOTEJID ] # runs in local queue
then
    job_id=$LSB_JOBID; array_id=$LSB_JOBINDEX
else # job has ben forwarded to remote queue
    job_id=$LSB_REMOTEJID; array_id=$LSB_REMOTEINDEX
fi

# 13 lambdas: 0 to 12
lambda_index=$(( ($array_id - 1) / 5 ))
# replicas from 1 to 5
replica=$(( $array_id - $lambda_index * 5 ))
echo "lambda $lambda_index"
echo "replica $replica"

# look up the lambda in a list
declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
# declare -a replicas=(1 2 3 4 5)
lambda=${lambdas[$lambda_index]}
echo "lambda $lambda"

cmd_mpinamd="mpiexec -n $NP namd2 +pemap 1-31 +commap 0"

# check the ligand first
cd lig
cd lambda_$lambda/rep$replica
# lig min
if test -f min.log && grep -q "WRITING VELOCITIES TO OUTPUT FILE" min.log; then
    echo "LIG MIN ALREADY DONE"
else
    echo "LIG MIN"
	${cmd_mpinamd} min.namd > min.log
fi
# eq1
if test -f eq_step1.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log; then
    echo "LIG EQ1 ALREADY DONE"
else
    echo "LIG EQ1"
	${cmd_mpinamd} eq_step1.namd > eq_step1.log
fi
# eq2
if test -f eq_step2.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step2.log; then
    echo "LIG EQ2 ALREADY DONE"
else
    echo "LIG EQ2"
	${cmd_mpinamd} eq_step2.namd > eq_step2.log
fi
# eq3
if test -f eq_step3.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step3.log; then
    echo "LIG EQ3 ALREADY DONE"
else
    echo "LIG EQ3"
	${cmd_mpinamd} eq_step3.namd > eq_step3.log
fi
# eq4
if test -f eq_step4.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step4.log; then
    echo "LIG EQ4 ALREADY DONE"
else
    echo "LIG EQ4"
	${cmd_mpinamd} eq_step4.namd > eq_step4.log
fi

# prod
if test -f prod.log && grep -q "WRITING VELOCITIES TO OUTPUT FILE AT STEP 3000000" prod.log; then
    echo "LIG PROD ALREADY DONE"
else
    echo "LIG PROD"
	${cmd_mpinamd} prod.namd > prod.log
fi
