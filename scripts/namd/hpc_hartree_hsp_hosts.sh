#!/bin/bash
#BSUB -J trl1l6
#BSUB -o std.%J.o
#BSUB -e std.%J.e
#BSUB -R "span[ptile=32]" ###### 32 cores per node
# this should be ptile * number of nodes,
#BSUB -n 2080       ##### for 64 cores, total nodes requested 2, 8320 is 128*65
#BSUB -q compbiomed
#BSUB -W 20:00
#BSUB -x

"""
This script is a simpler variation of the other one that has fewer dependancies.
"""

ROOT_WORK=$HCBASE/resp/tyk2_l1_l6
cd $ROOT_WORK
NP=32 # cores per simulation
HOSTS_PER_SIM=1 # numer of hosts per simulation
LIG_TIMEOUT=$(( 60 * 60 * 8 )) # time out of the ligand simulations
LIG_TIMEOUT=1

echo "Time start" `date`
env > last_used_env

#Load modules
source /etc/profile.d/modules.sh
module load namd-gcc/2.12

# generate all the simulations you want to run:
declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
declare -a replicas=(1 2 3 4 5)
# declare -a lambdas=(0.60) # testing
# declare -a replicas=(1 2) # testing
echo "Number lambdas: ${#lambdas[@]}"
echo "Number replicas: ${#replicas[@]}"

if [ -z $LSB_REMOTEJID ] # runs in local queue
then
    GRANTED_HOSTS=$LSB_HOSTS
else # job has ben forwarded to remote queue
    GRANTED_HOSTS=$LSB_MCPU_HOSTS
fi

env > envstd.%J

# print the number of threads
IFS=' ' read -r -a HOSTS_ARR <<< "${GRANTED_HOSTS}"
echo "Number of hyper thread cores: ${#HOSTS_ARR[@]}"

# get simulations
SIMS=()
for lambda in "${lambdas[@]}"; do
    for replica in "${replicas[@]}"; do
    	sys="lambda_$lambda/rep$replica"
    	SIMS+=($sys)
    done
done
SIM_NO=${#SIMS[@]}
echo "Number of simulations: $SIM_NO"
echo "Simulations: ${SIMS[@]}"


# get unique hosts
HOSTS=()
for host in ${GRANTED_HOSTS};
do
	# ignore the host if it's in the list already
	if [[ "${HOSTS[@]}" =~ "${host}" ]]; then
		continue
	fi
	# add a new host
    HOSTS+=($host)
done
echo "All hosts: ${HOSTS[@]}"


# schedule all simulations
for sim_no in $(seq 1 $SIM_NO); do
	# group hosts
	HOSTGROUP=" "
	for host in $(seq 1 $HOSTS_PER_SIM); do
		next_host=("${HOSTS[0]}") # get the first host
		HOSTS=("${HOSTS[@]:1}") # update the list
		# create the string for "mpirun -hosts VAR"
		if [ "$HOSTGROUP" != " " ]; then
        	HOSTGROUP="$HOSTGROUP,$next_host"
        else
        	HOSTGROUP="$next_host"
		fi
	done
	echo "Grouped Hosts for the next simulation: $HOSTGROUP"

	# get the next simulation
	SIM=("${SIMS[0]}")
	echo "Run the next $SIM"
	SIMS=("${SIMS[@]:1}") # remove it from the list

	# run the next simulation
	cmd_mpinamd="mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0"
	export I_MPI_PIN=yes
    (
        # example of a single line
        # echo hello && if [ $counter -eq 1 ]; then echo "true"; fi && echo done

	    # first is the ligand part with the timeout
	    ( cd $ROOT_WORK/lig/$SIM &&
	    timeout $LIG_TIMEOUT echo "Ligand simulation start" &&
	    if ! grep -q "WRITING VELOCITIES TO OUTPUT FILE" min.log ; then ${cmd_mpinamd} min.namd > min.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log ; then ${cmd_mpinamd} eq_step1.namd > eq_step1.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step2.log ; then ${cmd_mpinamd} eq_step2.namd > eq_step2.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step3.log ; then ${cmd_mpinamd} eq_step3.namd > eq_step3.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step4.log ; then ${cmd_mpinamd} eq_step4.namd > eq_step4.log ; fi &&
        if ! grep -q "WRITING VELOCITIES TO OUTPUT FILE AT STEP 3000000" prod.log ; then ${cmd_mpinamd} prod.namd > prod.log ; fi &&
        echo "Finished running lig/$SIM"
        ) ;
        (
        # then the complex part, no timeout
        cd $ROOT_WORK/complex/$SIM &&
        if ! grep -q "WRITING VELOCITIES TO OUTPUT FILE" min.log ; then ${cmd_mpinamd} min.namd > min.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log ; then ${cmd_mpinamd} eq_step1.namd > eq_step1.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log ; then ${cmd_mpinamd} eq_step2.namd > eq_step2.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log ; then ${cmd_mpinamd} eq_step3.namd > eq_step3.log ; fi &&
        if ! grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log ; then ${cmd_mpinamd} eq_step4.namd > eq_step4.log ; fi &&
        if ! grep -q "WRITING VELOCITIES TO OUTPUT FILE AT STEP 3000000" prod.log ; then ${cmd_mpinamd} prod.namd > prod.log ; fi &&
        echo "Finished running complex/$SIM"
        )
    ) &
done

echo "END: Sim that were not scheduled: ${SIMS[@]}"

wait
echo "Time End" `date`