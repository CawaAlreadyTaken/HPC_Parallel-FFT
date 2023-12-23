#!/bin/bash

#PBS -l select=1:ncpus=32:mem=2gb

#PBS -l walltime=0:01:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 32 HPC_Parallel-FFT/FFT_parallel/parallel_solver_1 HPC_Parallel-FFT/FFT_parallel/
