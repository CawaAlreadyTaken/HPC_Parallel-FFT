#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>

const int PRINTING_OUTPUT = 0;
const int PRINTING_TIME = 1;
MPI_Datatype mpi_send_tuple_type;

typedef struct {
    double real;
    double imag;
} complex;

typedef struct {
    complex value;
    int index;
} send_tuple;

int reverse(int num, int lg_n) {
    int res = 0;
    int i;
    for (i = 0; i < lg_n; i++) {
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

void multiply_transformed_polynomials(complex *a, int n0, complex *b, int n1) {
    assert(n0 == n1); // n0 and n1 must be equal: polynomials of the same length
	
    int i;
    for (i = 0; i < n0; i++) {
        a[i] = mul(a[i], b[i]);
    }
}

send_tuple * parallel_fft(complex *a, int n, int my_rank, int comm_sz, int lg_n, int invert) {
    // Calculating my start and end
    int my_start = my_rank * n / comm_sz;
    int my_end = (my_rank + 1) * n / comm_sz; // This is excluded
    int my_size = my_end - my_start;

    // Calculating lg_comm_sz
    int lg_comm_sz = 0;
    while ((1 << lg_comm_sz) < comm_sz) {
        lg_comm_sz++;
    }

    // This is the number of cycles for which we don't need to exchange data
    int no_exchange = lg_n - lg_comm_sz;

    // This is a counter for the current cycle
    int cycles = 1;

    // Allocate memory for the tuples we will send and receive
    send_tuple *to_send = malloc(my_size * sizeof(send_tuple));
    send_tuple *received = malloc(my_size * sizeof(send_tuple));

    int exchange_cycle;

    // Execute the first log(n)-log(comm_sz) cycles, for which we don't need to exchange data
    int len;
    for (len = 2; len <= n; len <<= 1) {
        double ang = 2*M_PI / len * (invert ? -1 : 1);
        complex wlen = complex_from_polar(1.0, ang);
        if (len > n/comm_sz) break;
        int i;
        for (i = my_start; i < my_end; i += len) {
            complex w = {1.0, 0.0};
            int j;
            for (j = 0; j < len / 2; j++) {
                complex u = a[i + j];
                complex v = mul(a[i + j + len / 2], w);
                a[i + j] = add(u, v);
                a[i + j + len / 2] = sub(u, v);
                w = mul(w, wlen);
            }
        }
        cycles++;
    }

    // Execute the last log(comm_sz) cycles, for which we need to exchange data
    for (len = 2*n/comm_sz; len <= n; len <<= 1) {
        exchange_cycle = cycles - no_exchange;
        double ang = 2*M_PI / len * (invert ? -1 : 1);
        complex wlen = complex_from_polar(1.0, ang);
        int send_index = 0;
        int number_of_iterations = 0;
        int i;
        for (i = 0; i < n; i += len) {
            complex w = {1.0, 0.0};
            int j;
            for (j = 0; j < len / 2; j++) {
                if (number_of_iterations >= my_end)
                    break;
                if (number_of_iterations >= my_start) {
                    complex u = a[i + j];
                    complex v = mul(a[i + j + len / 2], w);
                    a[i + j] = add(u, v);
                    a[i + j + len / 2] = sub(u, v);
                    w = mul(w, wlen);

                    // Fill the tuple array with the values we need to send
                    to_send[send_index].value = a[i+j];
                    to_send[send_index].index = i+j;
                    send_index++;
                    to_send[send_index].value = a[i+j + (len / 2)];
                    to_send[send_index].index = i+j + (len/2);
                    send_index++;
                }
                number_of_iterations+=2;
            }
        }

        // Distance between threads that will communicate with each other
        int distance = pow(2, exchange_cycle);

        // Last cycle we don't have anything to send to others
        if (len == n){
            free(received);
            // If needed apply the inverse transform
            if (invert) {
                int i;
                for (i = my_start; i < my_end; i++) {
                    a[i].real /= n;
                    a[i].imag /= n;
                }
            }
            return to_send;
        }
        if ((my_rank / distance) % 2 == 0){
            //printf("rank: %d, receiving from %d\n", my_rank, my_rank + distance);
            MPI_Recv(received, my_size, mpi_send_tuple_type, my_rank + distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //printf("rank: %d, received from %d\n", my_rank, my_rank + distance);
            //printf("rank: %d, sending to %d\n", my_rank, my_rank + distance);
            MPI_Send(to_send, my_size, mpi_send_tuple_type, my_rank + distance, 0, MPI_COMM_WORLD);
            //printf("rank: %d, sent to %d\n", my_rank, my_rank + distance);
        } else {
            //printf("rank: %d, sending to %d\n", my_rank, my_rank - distance);
            MPI_Send(to_send, my_size, mpi_send_tuple_type, my_rank - distance, 0, MPI_COMM_WORLD);
            //printf("rank: %d, sent to %d\n", my_rank, my_rank - distance);
            //printf("rank: %d, receiving from %d\n", my_rank, my_rank - distance);
            MPI_Recv(received, my_size, mpi_send_tuple_type, my_rank - distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //printf("rank: %d, received from %d\n", my_rank, my_rank - distance);
        }

        int x;
        for (x=0; x<my_size; x++){
            a[received[x].index] = received[x].value;
        }
        cycles++;
    }

    free(received);
    // If needed apply the inverse transform
    if (invert) {
	int i;
	for (i = my_start; i < my_end; i++) {
	    a[i].real /= n;
	    a[i].imag /= n;
	}
    }
    return NULL;
}

void gather_data(send_tuple * to_send, int my_size, int my_rank, complex * a, int n){
    if (my_rank == 0){
	if (to_send == NULL) return; // This happens when the comm_sz is 1. We already have all the data
        send_tuple* final_receive = malloc(sizeof(send_tuple) * n);
        MPI_Gather(to_send, my_size, mpi_send_tuple_type, final_receive, my_size, mpi_send_tuple_type, 0, MPI_COMM_WORLD);
        int x;
        for (x=0; x<n; x++){
            a[final_receive[x].index] = final_receive[x].value;
        }
    } else {
        MPI_Gather(to_send, my_size, mpi_send_tuple_type, NULL, my_size, mpi_send_tuple_type, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Timers declaration for measuring time
    clock_t start, end;

    // MPI new type definition
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    MPI_Aint offsets[2];

    offsets[0] = offsetof(send_tuple, value);
    offsets[1] = offsetof(send_tuple, index);

    MPI_Type_create_struct(2, blocklengths, offsets, types, &mpi_send_tuple_type);
    MPI_Type_commit(&mpi_send_tuple_type);

    if (my_rank == 0) {
        // We need to know the working directory
        if (argc != 2) {
            printf("Usage: %s <WORK_DIR>\n", argv[0]);
            printf("Example: %s ./\n", argv[0]);
            return 1;
        }

        // Opening file for writing time results
        const char *timings_file_name = "timing_parallel_solver_1.txt";
        int timings_file_length = strlen(argv[1]) + strlen(timings_file_name) + 1;
        char *full_timings_file = (char *)malloc(timings_file_length);
        strcpy(full_timings_file, argv[1]);
        strcat(full_timings_file, timings_file_name);
        FILE *timings_file = fopen(full_timings_file, "w");
        // Opening file for reading input
        const char *input_file_name = "../dataset/data/dataset_1_2.txt";
        int input_file_length = strlen(argv[1]) + strlen(input_file_name) + 1;
        char *full_input_file = (char *)malloc(input_file_length);
        strcpy(full_input_file, argv[1]);
        strcat(full_input_file, input_file_name);
        FILE *input_file = fopen(full_input_file, "r");
        if (timings_file == NULL) {
            fprintf(stderr, "Timings file was: %s\n", full_timings_file);
            perror("Error creating timings file");
            return 1;
        }
        if (input_file == NULL) {
            fprintf(stderr, "Input file was: %s\n", full_input_file);
            perror("Error opening input file");
            return 1;
        }
        free(full_timings_file);
        free(full_input_file);

        // Reading first input size
        int n0;
        fscanf(input_file, "%d", &n0);

        start = clock();

        // Allocating memory for first input array
        complex *a = malloc(n0 * sizeof(complex));

        // Reading first input array
        int i;
        for (i = 0; i < n0; i++) {
            fscanf(input_file, "%lf", &a[i].real);
            a[i].imag = 0;
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for allocating memory and reading first input: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Reordering the array
        int lg_n = 0;
        while ((1 << lg_n) < n0)
            lg_n++;

        // TODO: check data dependencies
        for (i = 0; i < n0; i++) {
            if (i < reverse(i, lg_n))
                swap(&a[i], &a[reverse(i, lg_n)]);
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for reordering the first input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Broadcast n0
        MPI_Bcast(&n0, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Broadcast the first array a
        MPI_Bcast(a, n0, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for broadcasting the first input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        int my_end_a = n0 / comm_sz;
        int my_size_a = my_end_a;

        start = clock();

        // Solve my data, retrieve the data to send in the last iteration
        send_tuple* data_to_send_a = parallel_fft(a, n0, my_rank, comm_sz, lg_n, 0);

        end = clock();
        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for calculating the first parallel fft: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Gather result in the root node
        gather_data(data_to_send_a, my_size_a, my_rank, a, n0);

        end = clock();
        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for first gathering data: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        // Reading second input size
        int n1;
        fscanf(input_file, "%d", &n1);

        start = clock();

        // Allocating memory for second array
        complex *b = malloc(n1 * sizeof(complex));

        // Reading input array
        for (i = 0; i < n1; i++) {
            fscanf(input_file, "%lf", &b[i].real);
            b[i].imag = 0;
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for allocating memory and reading second input: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }
        fclose(input_file);

        start = clock();

        // Reordering the second array
        lg_n = 0;
        while ((1 << lg_n) < n0)
            lg_n++;

        for (i = 0; i < n1; i++) {
            if (i < reverse(i, lg_n))
                swap(&b[i], &b[reverse(i, lg_n)]);
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for reordering the second input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Broadcast n1
        MPI_Bcast(&n1, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Broadcast the second array b
        MPI_Bcast(b, n1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for broadcasting the second input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        int my_end_b = n1 / comm_sz;
        int my_size_b = my_end_b;

        start = clock();

        // Solve my data, retrieve the data to send in the last iteration
        send_tuple* data_to_send_b = parallel_fft(b, n1, my_rank, comm_sz, lg_n, 0);

        end = clock();
        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for calculating the second parallel fft: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Gather result in the root node
        gather_data(data_to_send_b, my_size_b, my_rank, b, n1);

        end = clock();
        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for gathering second data: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Multiply polynomials
        multiply_transformed_polynomials(a, n0, b, n1);

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for multiplying polynomials: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        if (PRINTING_OUTPUT) {
            start = clock();

            printf("FFT result:\n");
            for (i = 0; i < n0; i++) {
                printf("(%f, %f)\n", a[i].real, a[i].imag);
            }

            end = clock();
            if (PRINTING_TIME) {
                fprintf(timings_file, "Time for printing output: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
            }
        }

        free(a);
        free(b);
        free(data_to_send_a);
        free(data_to_send_b);
    } else {
        int n0;
        MPI_Bcast(&n0, 1, MPI_INT, 0, MPI_COMM_WORLD);

        complex *a = malloc(n0 * sizeof(complex));
        MPI_Bcast(a, n0, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        int lg_n = 0;
        while ((1 << lg_n) < n0)
            lg_n++;

        int my_start_a = my_rank * n0 / comm_sz;
        int my_end_a = (my_rank + 1) * n0 / comm_sz;
        int my_size_a = my_end_a - my_start_a;

        send_tuple* data_to_send_a = parallel_fft(a, n0, my_rank, comm_sz, lg_n, 0);

	// Send data to the root node
        gather_data(data_to_send_a, my_size_a, my_rank, a, n0);

        int n1;
        MPI_Bcast(&n1, 1, MPI_INT, 0, MPI_COMM_WORLD);

        complex *b = malloc(n1 * sizeof(complex));
        MPI_Bcast(b, n1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        lg_n = 0;
        while ((1 << lg_n) < n1)
            lg_n++;

        int my_start_b = my_rank * n1 / comm_sz;
        int my_end_b = (my_rank + 1) * n1 / comm_sz;
        int my_size_b = my_end_b - my_start_b;

        send_tuple* data_to_send_b = parallel_fft(b, n1, my_rank, comm_sz, lg_n, 0);

	// Send data to the root node
        gather_data(data_to_send_b, my_size_b, my_rank, b, n1);

        free(a);
        free(b);
        free(data_to_send_a);
        free(data_to_send_b);
    }

    MPI_Type_free(&mpi_send_tuple_type);
    MPI_Finalize();
    return 0;
}
