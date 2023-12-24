#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
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

send_tuple * parallel_fft(complex *a, int n, int my_rank, int comm_sz, int lg_n, int invert) {
	// Calculating my start and end
	int my_start = n / comm_sz * my_rank;
	int my_end = n / comm_sz * (my_rank + 1); // This is excluded
	int my_size = my_end - my_start;

	// Calculating lg_comm_sz
	int lg_comm_sz = 0;
	while ((1 << lg_comm_sz) < comm_sz)
		lg_comm_sz++;

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
			printf("rank: %d, receiving from %d\n", my_rank, my_rank + distance);
			MPI_Recv(received, my_size, mpi_send_tuple_type, my_rank + distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank: %d, received from %d\n", my_rank, my_rank + distance);
			printf("rank: %d, sending to %d\n", my_rank, my_rank + distance);
			MPI_Send(to_send, my_size, mpi_send_tuple_type, my_rank + distance, 0, MPI_COMM_WORLD);
			printf("rank: %d, sent to %d\n", my_rank, my_rank + distance);
		} else {
			printf("rank: %d, sending to %d\n", my_rank, my_rank - distance);
			MPI_Send(to_send, my_size, mpi_send_tuple_type, my_rank - distance, 0, MPI_COMM_WORLD);
			printf("rank: %d, sent to %d\n", my_rank, my_rank - distance);
			printf("rank: %d, receiving from %d\n", my_rank, my_rank - distance);
			MPI_Recv(received, my_size, mpi_send_tuple_type, my_rank - distance, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank: %d, received from %d\n", my_rank, my_rank - distance);
		}

		int x;
		for (x=0; x<my_size; x++){
			a[received[x].index] = received[x].value;
		}
		cycles++;
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

int main(int argc, char* argv[]) {
	int comm_sz;
	int my_rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Timers declaration for measuring time
	clock_t start, end;

	// MPI new type definition
	int blocklengths[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_DOUBLE_COMPLEX, MPI_INT};
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
		const char *timings_file_name = "timing_parallel_solver_0.txt";
		int timings_file_length = strlen(argv[1]) + strlen(timings_file_name) + 1;
		char *full_timings_file = (char *)malloc(timings_file_length);
		strcpy(full_timings_file, argv[1]);
		strcat(full_timings_file, timings_file_name);
		FILE *timings_file = fopen(full_timings_file, "w");
		// Opening file for reading input
		const char *input_file_name = "../dataset/data/dataset_0_3.txt";
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

		// Reading input size
		int n;
		fscanf(input_file, "%d", &n);

		start = clock();

		// Allocating memory for input array
		complex *input_a = malloc(n * sizeof(complex));
		complex *a = malloc(n * sizeof(complex));

		// Reading input array
		int i;
		for (i = 0; i < n; i++) {
			fscanf(input_file, "%lf", &input_a[i].real);
			input_a[i].imag = 0;
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

		// TODO: check data dependencies
		for (i = 0; i < n; i++) {
			int rev = reverse(i, lg_n);
			if (i < rev)
				swap(&input_a[i], &input_a[rev]);
		}

		end = clock();

		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for reordering the input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		start = clock();

		// Broadcast n
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

		int my_end = n / comm_sz;
		int my_size = my_end;

		// Scatter data
		MPI_Scatter(input_a, my_size, mpi_send_tuple_type, a, my_size, mpi_send_tuple_type, 0, MPI_COMM_WORLD);

		end = clock();

		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for broadcasting the input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		start = clock();

		// Solve my data
		send_tuple* data_to_send = parallel_fft(a, n, my_rank, comm_sz, lg_n, 0);

		end = clock();
		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for calculating the parallel fft: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		start = clock();

		// Gather result in the root node
		gather_data(data_to_send, my_size, my_rank, a, n);

		end = clock();
		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for gathering data: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		if (PRINTING_OUTPUT) {
			start = clock();

			// Print the result
			for (i=0; i<n; i++){
				printf("%lf ", a[i].real);
			}
			printf("\n");

			end = clock();
			if (PRINTING_TIME) {
				fprintf(timings_file, "Time for printing result: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
			}
		}

		fclose(timings_file);

		// Free memory
		free(a);
	} else {
		int n;
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

		int my_start = n / comm_sz * my_rank;
		int my_end = n / comm_sz * (my_rank + 1);
		int my_size = my_end - my_start;

		complex *a = malloc(n * sizeof(complex));
		MPI_Scatter(NULL, my_size, mpi_send_tuple_type, &a[my_start], my_size, mpi_send_tuple_type, 0, MPI_COMM_WORLD);

		int lg_n = 0;
		while ((1 << lg_n) < n)
			lg_n++;

		send_tuple* data_to_send = parallel_fft(a, n, my_rank, comm_sz, lg_n, 0);

		gather_data(data_to_send, my_size, my_rank, a, n);

		free(a);
	}

	MPI_Type_free(&mpi_send_tuple_type);
	MPI_Finalize();
	return 0;
}
