#!/bin/bash
# define an array job
#BSUB -J mr5_32_38[1-65]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 96
##BSUB -q compbiomed
#BSUB -q scafellpikeSKL # can be used a litte, 10x4, try 65x1 (does not work). Seems to allow only for 10 at a time.
#BSUB -W 35:00
#BSUB -x

# schedule 5 simulaitons for each instead,
# use the variables to determine which host to use, like in the other scripts
# you can try running 10 x 5 to get 50 sims done in one go, which would be handy

SYSTEM=$HCBASE/resp_ff/mcl1_l32_l38_5th
NP=96

#Load modules
source /etc/profile.d/modules.sh
module load namd-gcc/2.12

# look up the index for this job instance
if [ -z $LSB_REMOTEJID ] # runs in local queue
then
    job_id=$LSB_JOBID; array_id=$LSB_JOBINDEX
    GRANTED_HOSTS=$LSB_HOSTS
else # job has ben forwarded to remote queue
    job_id=$LSB_REMOTEJID; array_id=$LSB_REMOTEINDEX
    GRANTED_HOSTS=$LSB_MCPU_HOSTS
fi

# print the number of threads
IFS=' ' read -r -a HOSTS_ARR <<< "${GRANTED_HOSTS}"
echo "Number of hyper thread cores: ${#HOSTS_ARR[@]}"

# get unique hosts
HOSTS=()
for host in ${GRANTED_HOSTS};
do
    # ignore a host that is a number
    re='^[0-9]+$'
    if [[ $host =~ $re ]] ; then
       continue
    fi

	# ignore the host if it's in the list already
	if [[ "${HOSTS[@]}" =~ "${host}" ]]; then
		continue
	fi
	# add a new host
    HOSTS+=($host)
    echo "Adding new host: $host"
done
echo "All hosts: ${HOSTS[@]}"
# fixme - add checks that hosts have been extracted, throw an error otherwise

# so for each host you need keep a counter and continue to run simulation?
# start from replica 1
replica=1
for host in ${HOSTS};
do
    echo "Next host is $host"

    # 13 lambdas: 0 to 12
    lambda_index=$(( ($array_id - 1) / 5 ))
    echo "lambda $lambda_index"
    echo "replica $replica"
    env > out%J_env_$array_id

    # look up the lambda in a list
    declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
    # declare -a replicas=(1 2 3 4 5)
    lambda=${lambdas[$lambda_index]}
    echo "lambda $lambda"

    cmd_mpinamd="mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $host namd2 +pemap 1-31 +commap 0"
    (
        # check the ligand first
        echo "Running the first system"
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
    ) &

    # switch to the next replica
    replica=$(( $replica + 1 ))
    echo "Increasing replica. Now value: $replica"
done