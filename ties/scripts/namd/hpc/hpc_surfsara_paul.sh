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
GROMACS/2019.3-intel-2018b-CUDA-10.0.130

# energy minimisation
mpirun -np $SLURM_CPUS_ON_NODE namd2 min.namd > min.log
# equilibriate with different constraints
