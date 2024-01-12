# HPC_Parallel-FFT

## Serial Solution
The serial solution is in `serial_solver/`.  

`serial_solver/serial_solver_0.c` solves the task of just calculating the DFT transform of one polynomial.  

`serial_solver/serial_solver_1.c` solves the task of multiplying two polynomials using the FFT.

## Parallel Solution
The parallel solution is in `FFT_parallel`.  

Again, `FFT_parallel/parallel_solver_0.c` solves the first task, and `FFT_parallel/parallel_solver_1.c` solves the second task.  

`FFT_parallel/pack_excl_timings_0.txt` and `FFT_parallel/pack_excl_timings_1.txt` hold the raw timings that are published in the report.  
