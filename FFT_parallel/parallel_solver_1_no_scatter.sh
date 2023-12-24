#!/bin/bash

#PBS -l select=1:ncpus=8:mem=16gb

#PBS -l walltime=0:10:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 8 HPC_Parallel-FFT/FFT_parallel/parallel_solver_1_no_scatter HPC_Parallel-FFT/FFT_parallel/
