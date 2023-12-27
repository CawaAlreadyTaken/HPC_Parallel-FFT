#!/bin/bash

#PBS -l select=1:ncpus=16:mem=8gb -l place=pack:excl

#PBS -l walltime=1:50:00

#PBS -q short_cpuQ

#PBS -e easd
#PBS -o oasd

module load mpich-3.2
mpirun.actual -n 16 HPC_Parallel-FFT/FFT_parallel/parallel_solver_1 HPC_Parallel-FFT/FFT_parallel/
