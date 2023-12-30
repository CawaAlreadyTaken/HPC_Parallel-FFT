#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(){
    int comm_sz;
	int my_rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int i;
    if (my_rank == 0)
    {
        #pragma omp for
        for (i=0; i<5; i++){
            printf("thread id: %d \n", omp_get_thread_num()); //print thread id to check if oprnmp works
        }
    }else{
        printf("I'm not rank 0 \n");
    }
    return 0;
}