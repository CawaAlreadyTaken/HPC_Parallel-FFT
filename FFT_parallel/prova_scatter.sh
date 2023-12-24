#!/bin/bash

#PBS -l select=4:ncpus=2:mem=2gb

#PBS -l walltime=0:02:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual --verbose -n 8 HPC_Parallel-FFT/FFT_parallel/prova_scatter HPC_Parallel-FFT/FFT_parallel/
