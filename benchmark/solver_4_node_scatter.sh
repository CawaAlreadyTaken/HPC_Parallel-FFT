#!/bin/bash

#PBS -l select=4:ncpus=1:mem=2gb place=scatter

#PBS -l walltime=0:02:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 32 HPC_Parallel-FFT/FFT_parallel/parallel_solver_0 HPC_Parallel-FFT/FFT_parallel/
