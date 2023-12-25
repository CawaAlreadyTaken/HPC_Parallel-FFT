#!/bin/bash

#PBS -l select=4:ncpus=16:mem=2gb

#PBS -l walltime=0:06:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 64 HPC_Parallel-FFT/FFT_parallel/parallel_solver_1 HPC_Parallel-FFT/FFT_parallel/
