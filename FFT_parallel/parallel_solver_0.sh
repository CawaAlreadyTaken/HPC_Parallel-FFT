#!/bin/bash

#PBS -l select=1:ncpus=16:mem=4gb -l place=pack:excl

#PBS -l walltime=0:05:00

#PBS -q short_cpuQ

#PBS -e easd
#PBS -o oasd

module load mpich-3.2
mpirun.actual -n 8 HPC_Parallel-FFT/FFT_parallel/parallel_solver_0 HPC_Parallel-FFT/FFT_parallel/
