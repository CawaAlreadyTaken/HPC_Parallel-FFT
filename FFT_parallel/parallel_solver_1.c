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

int min(int a, int b) {
    return (a < b ? a:b);
}

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
    int minimum = min(n0, n1);
    // We can just continue until minimum, the others will be 0
    int i;
    for (i = 0; i < minimum; i++) {
        a[i] = mul(a[i], b[i]);
    }
}

send_tuple * parallel_fft(complex *a, int n, int my_rank, int comm_sz, int lg_n, int invert) {
    // Calculating my start and end
    int my_start = n / comm_sz * my_rank;
    int my_end = n / comm_sz * (my_rank + 1); // This is excluded
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
            MPI_Recv(received, my_size, mpi_send_tuple_type, my_rank + distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(to_send, my_size, mpi_send_tuple_type, my_rank + distance, 0, MPI_COMM_WORLD);
        } else {
            MPI_Send(to_send, my_size, mpi_send_tuple_type, my_rank - distance, 0, MPI_COMM_WORLD);
            MPI_Recv(received, my_size, mpi_send_tuple_type, my_rank - distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

void custom_scatter(int my_rank, int comm_sz, int n, complex * a) {
	int size = n / comm_sz;
	if (my_rank == 0) {
		int i;
		for (i = 1; i < comm_sz; i++) {
			MPI_Send(a + i * size, size, MPI_DOUBLE_COMPLEX, i, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(a, size, MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

void Build_mpi_type(
	complex* 	value,
	int*		index
) {
	int array_of_blocklengths[2] = {1, 1};
	MPI_Datatype array_of_types[2] = {MPI_DOUBLE_COMPLEX, MPI_INT};
	MPI_Aint value_addr, index_addr;
	MPI_Aint array_of_displacements[2] = {0};
	MPI_Get_address(value, &value_addr);
	MPI_Get_address(index, &index_addr);
	array_of_displacements[1] = index_addr-value_addr;
	MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements, array_of_types, &mpi_send_tuple_type);
	MPI_Type_commit(&mpi_send_tuple_type);
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
    send_tuple x;
    Build_mpi_type(&x.value, &x.index);

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
        const char *input_file_name = "../dataset/data/dataset_1_4.txt";
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
	float t0;
	float t1;
	float t2;
	float t3;
	float t4;

        start = clock();
        // Reading first input size
        int n0;
        fscanf(input_file, "%d", &n0);


        // Allocating memory for first input array
        complex *a = malloc(2 * n0 * sizeof(complex));

        // Reading first input array
        int i;
        for (i = 0; i < n0; i++) {
            fscanf(input_file, "%lf", &a[i].real);
            a[i].imag = 0;
        }
	for (i = n0; i < 2*n0; i++) {
	    a[i].real = 0;
	    a[i].imag = 0;
	}

        end = clock();

        t0 = end-start;

        start = clock();

        // Reordering the array
        int lg_n = 0;
        while ((1 << lg_n) < 2*n0)
            lg_n++;

        // TODO: check data dependencies
        for (i = 0; i < 2*n0; i++) {
	    int rev = reverse(i, lg_n);
            if (i < rev)
                swap(&a[i], &a[rev]);
        }

        end = clock();

        t1 = end-start;

        start = clock();

        // Broadcast n0
        MPI_Bcast(&n0, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int my_end_a = 2*n0 / comm_sz;
        int my_size_a = my_end_a;

        // Scatter the first array a
	custom_scatter(my_rank, comm_sz, 2*n0, a);

        end = clock();

        t2 = end-start;

        start = clock();

        // Solve my data, retrieve the data to send in the last iteration
        send_tuple* data_to_send_a = parallel_fft(a, 2*n0, my_rank, comm_sz, lg_n, 0);

        end = clock();

	t3 = end-start;

        start = clock();

        // Gather result in the root node
        gather_data(data_to_send_a, my_size_a, my_rank, a, 2*n0);

        end = clock();
        
	t4 = end-start;

        start = clock();
        // Reading second input size
        int n1;
        fscanf(input_file, "%d", &n1);

        // Allocating memory for second array
        complex *b = malloc(2 * n1 * sizeof(complex));

        // Reading input array
        for (i = 0; i < n1; i++) {
            fscanf(input_file, "%lf", &b[i].real);
            b[i].imag = 0;
        }
	for (i = n1; i < 2*n1; i++) {
	    b[i].real = 0;
	    b[i].imag = 0;
	}

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for allocating memory and reading input: %f seconds\n", (double)(t0 + end - start) / CLOCKS_PER_SEC);
        }
        fclose(input_file);

        start = clock();

        // Reordering the second array
        lg_n = 0;
        while ((1 << lg_n) < 2*n1)
            lg_n++;

        for (i = 0; i < 2*n1; i++) {
	    int rev = reverse(i, lg_n);
            if (i < rev)
                swap(&b[i], &b[rev]);
        }

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for reordering the second input array: %f seconds\n", (double)(t1 + end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Broadcast n1
        MPI_Bcast(&n1, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int my_end_b = 2*n1 / comm_sz;
        int my_size_b = my_end_b;

        // Scatter the second array b
	custom_scatter(my_rank, comm_sz, 2*n1, b);

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for scattering the second input array: %f seconds\n", (double)(t2 + end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Solve my data, retrieve the data to send in the last iteration
        send_tuple* data_to_send_b = parallel_fft(b, 2*n1, my_rank, comm_sz, lg_n, 0);

        end = clock();
        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for calculating the second parallel fft: %f seconds\n", (double)(t3 + end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Gather result in the root node
        gather_data(data_to_send_b, my_size_b, my_rank, b, 2*n1);

        end = clock();
        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for gathering second data: %f seconds\n", (double)(t4 + end - start) / CLOCKS_PER_SEC);
        }

        start = clock();

        // Multiply polynomials
        multiply_transformed_polynomials(a, 2*n0, b, 2*n1);

        end = clock();

        if (PRINTING_TIME) {
            fprintf(timings_file, "Time for multiplying polynomials: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }

        if (PRINTING_OUTPUT) {
            start = clock();

            printf("FFT result:\n");
	    int minimum = min(2*n0, 2*n1);
            for (i = 0; i < minimum; i++) {
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

        int my_start_a = 2*n0 / comm_sz * my_rank;
        int my_end_a = 2*n0 / comm_sz * (my_rank + 1);
        int my_size_a = my_end_a - my_start_a;

        complex *a = malloc(2 * n0 * sizeof(complex));
	custom_scatter(my_rank, comm_sz, 2*n0, a);

        int lg_n = 0;
        while ((1 << lg_n) < 2*n0)
            lg_n++;

        send_tuple* data_to_send_a = parallel_fft(a, 2*n0, my_rank, comm_sz, lg_n, 0);

	// Send data to the root node
        gather_data(data_to_send_a, my_size_a, my_rank, a, 2*n0);

        int n1;
        MPI_Bcast(&n1, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int my_start_b = 2*n1 / comm_sz * my_rank;
        int my_end_b = 2*n1 / comm_sz * (my_rank + 1);
        int my_size_b = my_end_b - my_start_b;

        complex *b = malloc(2*n1 * sizeof(complex));
	custom_scatter(my_rank, comm_sz, 2*n1, b);

        lg_n = 0;
        while ((1 << lg_n) < 2*n1)
            lg_n++;

        send_tuple* data_to_send_b = parallel_fft(b, 2*n1, my_rank, comm_sz, lg_n, 0);

	// Send data to the root node
        gather_data(data_to_send_b, my_size_b, my_rank, b, 2*n1);

        free(a);
        free(b);
        free(data_to_send_a);
        free(data_to_send_b);
    }

    MPI_Type_free(&mpi_send_tuple_type);
    MPI_Finalize();
    return 0;
}
