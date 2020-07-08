#!/bin/bash
# define an array job
A='test'
#BSUB -J thr_7_3_knl%A[1-65]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=8]"
#BSUB -n 128
#BSUB -q scafellpikeKNL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
#BSUB -W 25:00
#BSUB -x

SYSTEM=$HCBASE/bcc_ff/thrombin/l7_l3

#Load modules
source /etc/profile.d/modules.sh
module load namd-knl/2.13

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
env > out%J_env_$array_id

# look up the lambda in a list
declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
# declare -a replicas=(1 2 3 4 5)
lambda=${lambdas[$lambda_index]}
echo "lambda $lambda"

# count how many processors are allocated
NP=0
for TOKEN in $LSB_HOSTS
do
   ((NP++))
done

echo "NP is $NP"

# In theory we are using 16 nodes. However, they do not seem to be physical nodes.
# even though the env shows that there is 16 different names for hosts.
cmd_mpinamd="charmrun ++verbose ++n $NP ++mpiexec-no-n ++remote-shell mpiexec \
    /lustre/scafellpike/local/apps/intel/namd/2.13-knl/bin/namd2 +setcpuaffinity \
    ++ppn 7 +pemap 0-127:8.7 +commap 7-127:8 "

# check the ligand first
cd $SYSTEM/lig
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

# check the ligand first
cd $SYSTEM/complex
cd lambda_$lambda/rep$replica
# lig min
if test -f min.log && grep -q "WRITING VELOCITIES TO OUTPUT FILE" min.log; then
    echo "COMPLEX MIN ALREADY DONE"
else
    echo "COMPLEX MIN"
	${cmd_mpinamd} min.namd > min.log
fi
# eq1
if test -f eq_step1.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log; then
    echo "COMPLEX EQ1 ALREADY DONE"
else
    echo "COMPLEX EQ1"
	${cmd_mpinamd} eq_step1.namd > eq_step1.log
fi
# eq2
if test -f eq_step2.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step2.log; then
    echo "COMPLEX EQ2 ALREADY DONE"
else
    echo "COMPLEX EQ2"
	${cmd_mpinamd} eq_step2.namd > eq_step2.log
fi
# eq3
if test -f eq_step3.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step3.log; then
    echo "COMPLEX EQ3 ALREADY DONE"
else
    echo "COMPLEX EQ3"
	${cmd_mpinamd} eq_step3.namd > eq_step3.log
fi
# eq4
if test -f eq_step4.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step4.log; then
    echo "COMPLEX EQ4 ALREADY DONE"
else
    echo "COMPLEX EQ4"
	${cmd_mpinamd} eq_step4.namd > eq_step4.log
fi
# prod
if test -f prod.log && grep -q "WRITING VELOCITIES TO OUTPUT FILE AT STEP 3000000" prod.log; then
    echo "COMPLEX PROD ALREADY DONE"
else
    echo "COMPLEX PROD"
	${cmd_mpinamd} prod.namd > prod.log
fi