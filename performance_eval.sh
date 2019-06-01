#!/bin/bash -l
# The -l above is required to get the full environment with modules

# The name of the script is myjob
#SBATCH -J performance_eval

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 1:00:00
#SBATCH -A edu19.DD2356

# Number of nodes
#SBATCH --nodes=8

# Number of MPI processes per node
#SBATCH --ntasks-per-node=32

# Number of cores hosting OpenMP threads
#SBATCH -c 8

#SBATCH -e error_file.e

# TODO add threads
# export OMP_NUM_THREADS=8 


# Strong scalling, 20 GB

aprun -n 32 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_32p_20GB.txt

aprun -n 64 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_64p_20GB.txt

aprun -n 128 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_128p_20GB.txt

aprun -n 256 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt > ./results/results_256p_20GB.txt


# Strong scalling, 160 GB

aprun -n 32 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > ./results/results_32p_160GB.txt

aprun -n 64 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > ./results/results_64p_160GB.txt

aprun -n 128 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > ./results/results_128p_160GB.txt

aprun -n 256 -N 32 ./bin/project.out -r 5 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_160GB.txt > ./results/results_256p_160GB.txt