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
    int i;
    #pragma omp master
    {
        #pragma omp for num_threads(5)
        for (i=0; i<5; i++){
            printf("thread id: %d \n", omp_get_thread_num()); //print thread id to check if oprnmp works
        }
    }
    return 0;
}