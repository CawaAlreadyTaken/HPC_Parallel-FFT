#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

const double PI = acos(-1);
const int PRINTING_OUTPUT = 0;
const int PRINTING_TIME = 1;

typedef struct {
    double real;
    double imag;
} complex;

int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

void swap(complex *a, complex *b) {
    complex temp = *a;
    *a = *b;
    *b = temp;
}

complex complex_from_polar(double r, double theta_rad) {
    complex result;
    result.real = r * cos(theta_rad);
    result.imag = r * sin(theta_rad);
    return result;
}

complex add(complex a, complex b) {
    complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

complex sub(complex a, complex b) {
    complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

complex mul(complex a, complex b) {
    complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

void partial_fft(complex *a, int n, int my_rank, int comm_sz, int lg_n, int invert) {
    // Calculating my start and end
    int my_start = my_rank * n / comm_sz;
    int my_end = (my_rank + 1) * n / comm_sz; // This is excluded
    int my_size = my_end - my_start;

    // Calculating lg_comm_sz
    int lg_comm_sz = 0;
    while ((1 << lg_comm_sz) < comm_sz)
        lg_comm_sz++;

    // This is the number of cycles for which we don't need to exchange data
    int number_of_cycles = lg_n - lg_comm_sz;
    printf("number_of_cycles = %d\n", number_of_cycles);

    // Partial fft for the cycles for which we don't need to exchange data
    for (int len = 2; len <= n && number_of_cycles--; len <<= 1) {
        printf("partial number_of_cycles = %d\n", number_of_cycles);
        double ang = 2*PI / len * (invert ? -1 : 1);
        complex wlen = complex_from_polar(1.0, ang);
        for (int i = my_start; i < my_end; i += len) {
            complex w = {1.0, 0.0};
            for (int j = 0; j < len / 2; j++) {
                complex u = a[i + j];
                complex v = mul(a[i + j + len / 2], w);
                a[i + j] = add(u, v);
                a[i + j + len / 2] = sub(u, v);
                w = mul(w, wlen);
            }
        }
    }

    /* Only do this ad the end of fft, not when partial
    if (invert) {
        for (int i = my_start; i < my_end; i++) {
            a[i].real /= n;
            a[i].imag /= n;
        }
    }
    */

}

int main(void) {
    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Timers declaration for measuring time
    clock_t start, end;
    
    if (my_rank == 0) {
        // Opening file for writing time results
        FILE *timings_file = fopen("timing_parallel_solver_0.txt", "w");
        // Opening file for reading input
        FILE *input_file = fopen("input.txt", "r");
        if (timings_file == NULL) {
            perror("Error creating timing_parallel_solver_0.txt file.");
            return 1;
        }
        if (input_file == NULL) {
            perror("Error opening input.txt file.");
            return 1;
        }

        // Reading input size
        int n;
        fscanf(input_file, "%d", &n);

        start = clock();

        // Allocating memory for input array
        complex *a = malloc(n * sizeof(complex));

        // Reading input array
        for (int i = 0; i < n; i++) {
            fscanf(input_file, "%lf", &a[i].real);
            a[i].imag = 0;
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for allocating memory and reading input: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }
        fclose(input_file);

        start = clock();

        // Reordering the array
        int lg_n = 0;
        while ((1 << lg_n) < n)
            lg_n++;

        for (int i = 0; i < n; i++) {
            if (i < reverse(i, lg_n))
                swap(&a[i], &a[reverse(i, lg_n)]);
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for reordering the input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        //Broadcast n
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Broadcast the array a
        MPI_Bcast(a, n, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        //Solve my data
        partial_fft(a, n, my_rank, comm_sz, lg_n, 0);


        //Print my result


        fclose(timings_file);

        //Free memory
        free(a);
    } else {
        // Receiving n
        int n;
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Calculating lg_n
        int lg_n = 0;
        while ((1 << lg_n) < n)
            lg_n++;

        // Receiving the array a
        complex *a = malloc(n * sizeof(complex));
        MPI_Bcast(a, n, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        //Solve my data
        partial_fft(a, n, my_rank, comm_sz, lg_n, 0);


        //Print my result


        //Free memory
        free(a);
    }

    MPI_Finalize();
    return 0;
}
