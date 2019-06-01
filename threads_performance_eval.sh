#!/bin/bash

# The name of the script is myjob
#SBATCH -J threads_performance_eval

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 1:00:00
#SBATCH -A edu19.DD2356

# Number of nodes
#SBATCH --nodes=4

# Number of MPI processes per node
#SBATCH --ntasks-per-node=32

#SBATCH -e error_file.e

export OMP_NUM_THREADS=8
aprun -n 16 -N 16 -d 8 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_16p_8t_20GB.txt

export OMP_NUM_THREADS=4
aprun -n 32 -N 32 -d 4 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_32p_4t_20GB.txt

export OMP_NUM_THREADS=2
aprun -n 64 -N 32 -d 2 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_64p_2t_20GB.txt

export OMP_NUM_THREADS=1
aprun -n 128 -N 32 -d 1 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_128p_1t_20GB.txt

