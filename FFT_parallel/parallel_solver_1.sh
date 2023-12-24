#!/bin/bash

#PBS -l select=1:ncpus=1:mem=2gb

#PBS -l walltime=0:03:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 1 HPC_Parallel-FFT/FFT_parallel/parallel_solver_1 HPC_Parallel-FFT/FFT_parallel/
