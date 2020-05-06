#!/bin/bash
#BSUB -J t2rl5l16
#BSUB -o std.%J.o
#BSUB -e std.%J.e
#BSUB -R "span[ptile=32]" ###### 32 cores per node
# 64 cpu * 65 nodes
# this should be ptile * number of nodes,
#BSUB -n 2080       ##### for 64 cores, total nodes requested 2, 8320 is 128*65
#BSUB -q compbiomed
#BSUB -W 55:00
#BSUB -x

ROOT_WORK=$HCBASE/resp/tyk2_l5_l16
cd $ROOT_WORK
NP=32 # cores per simulation
HOSTS_PER_SIM=1 # numer of hosts per simulation
LIG_TIMEOUT=$(( 60 * 60 * 15 )) # time out of the ligand simulations


#Load modules
source /etc/profile.d/modules.sh
module load namd-gcc/2.12

# generate all the simulations you want to run:
declare -a lambdas=(0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00)
declare -a replicas=(1 2 3 4 5)
# declare -a lambdas=(0.50) # testing
# declare -a replicas=(1 2) # testing
echo "Number lambdas: ${#lambdas[@]}"
echo "Number replicas: ${#replicas[@]}"

# print the number of threads
IFS=' ' read -r -a HOSTS_ARR <<< "${LSB_HOSTS}"
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
#LSB_HOSTS="sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6b135.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx sqg6e4.bullx"
HOSTS=()
for host in ${LSB_HOSTS};
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
	export I_MPI_PIN=yes
    (

	    # first is the ligand part with the timeout
	    ( cd $ROOT_WORK/lig/$SIM &&
	    timeout $LIG_TIMEOUT mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 min.namd > min.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step1.namd > eq_step1.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step2.namd > eq_step2.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step3.namd > eq_step3.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step4.namd > eq_step4.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 prod.namd > prod.log &&
        echo "Finished running lig/$SIM"
        ) ;
        (
        # then the complex part, no timeout
        cd $ROOT_WORK/complex/$SIM &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 min.namd > min.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step1.namd > eq_step1.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step2.namd > eq_step2.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step3.namd > eq_step3.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 eq_step4.namd > eq_step4.log &&
        mpiexec -genv I_MPI_PIN_DOMAIN=auto:compact -n $NP -hosts $HOSTGROUP namd2 +pemap 1-31 +commap 0 prod.namd > prod.log &&
        echo "Finished running complex/$SIM"
        )
    ) &
done

echo "END: Sim that were not scheduled: ${SIMS[@]}"

wait