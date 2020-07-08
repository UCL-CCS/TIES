#!/bin/sh
# just one normal node
#SBATCH -N 1        # request 1 nodes
##SBATCH -p short    # request partition 'short', see below
##SBATCH -n 40       # request 40 processes (each runs on 1 core), mostly to run MPI programs on
                    # SLURM will compute the number of nodes needed
##SBATCH -n 16 -c 4  # request 16*4 cores, an MPI program will start
                    # 16 processes, each process can spawn 4 OpenMP
		    # threads. (the environment variable OMP_NUM_THREADS
		    # will be set to 4)
##SBATCH -N 6 --ntasks-per-node=4
                    # An MPI program will start 4 processes on 6 nodes each
		    # 24 processes in total
##SBATCH -t 2:00:00  # The job can take at most 2 wall-clock hours.
##SBATCH -t 10       # 10 minutes
##SBATCH -t 10:20    # 10 minutes plus 20 seconds
#SBATCH -t 15:20:30 #  10 hours plus 20 minutes plus 30 seconds

module load 2019
module load NAMD/2.13-foss-2018b-mpi
#module load pre2019
#module load NAMD/2.12-foss-2017b-mpi

# energy minimisation
mpirun -np $SLURM_CPUS_ON_NODE namd2 min.namd > min.log
# equilibriate with different constraints
mpirun -np $SLURM_CPUS_ON_NODE namd2 eq_step1.namd > eq_step1.log
mpirun -np $SLURM_CPUS_ON_NODE namd2 eq_step2.namd > eq_step2.log
mpirun -np $SLURM_CPUS_ON_NODE namd2 eq_step3.namd > eq_step3.log
mpirun -np $SLURM_CPUS_ON_NODE namd2 eq_step4.namd > eq_step4.log
# run the production
mpirun -np $SLURM_CPUS_ON_NODE namd2 prod.namd > prod.log