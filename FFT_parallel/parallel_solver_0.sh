#!/bin/bash

#PBS -l select=8:ncpus=4:mem=33gb

#PBS -l walltime=0:10:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual --verbose -n 32 HPC_Parallel-FFT/FFT_parallel/parallel_solver_0 HPC_Parallel-FFT/FFT_parallel/
