#!/bin/bash

#PBS -l select=2:ncpus=4:mem=2gb

#PBS -l walltime=0:01:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 8 HPC_Parallel-FFT/FFT_parallel/parallel_solver_1 HPC_Parallel-FFT/FFT_parallel/
