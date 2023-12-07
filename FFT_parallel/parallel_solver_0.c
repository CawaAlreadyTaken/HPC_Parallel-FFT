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
	double real;
	double imag;
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

void parallel_fft(complex *a, int n, int my_rank, int comm_sz, int lg_n, int invert) {
	// Calculating my start and end
	int my_start = my_rank * n / comm_sz;
	int my_end = (my_rank + 1) * n / comm_sz; // This is excluded
	int my_size = my_end - my_start;

	// Calculating lg_comm_sz
	int lg_comm_sz = 0;
	while ((1 << lg_comm_sz) < comm_sz)
		lg_comm_sz++;

	// This is the number of cycles for which we don't need to exchange data
	int no_exchange = lg_n - lg_comm_sz;

	//This is a counter for th cycle we are inside
	int cycles = 1;

	// Parallel fft
	int len;
	//send_tuple *to_send = malloc(my_size * sizeof(send_tuple)); //data to send to other thred --> data modified in this cycle !!CHECK LENGHT!!
	send_tuple to_send[my_size];
	//send_tuple *received = malloc(my_size * sizeof(send_tuple));
	send_tuple received[my_size];

	int exchange_cycle = 0;

	for (len = 2; len <= n; len <<= 1) {
		exchange_cycle = cycles - no_exchange;
		int send_index = 0;
		double ang = 2*M_PI / len * (invert ? -1 : 1);
		complex wlen = complex_from_polar(1.0, ang);
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
				/*
				if (exchange_cycle >=0) {
					to_send[send_index].real = a[i+j].real;
					to_send[send_index].imag = a[i+j].imag;
					to_send[send_index].index = i+j;
					send_index++;
					to_send[send_index].real = a[i+j + (len / 2)].real;
					to_send[send_index].imag = a[i+j + (len / 2)].imag;
					to_send[send_index].index = i+j + (len/2);
					send_index++;
				}
				*/
			}
		}

		/*
		if (exchange_cycle >= 0){
			int distance = pow(2, exchange_cycle); //distance between thread that will communicate with each other

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
			printf("rank: %d reached barrier at cycle %d\n", my_rank, exchange_cycle);
			MPI_Barrier(MPI_COMM_WORLD);
			printf("rank: %d passed barrier at cycle %d\n", my_rank, exchange_cycle);

			int x;
			for (x=0; x<my_size; x++){
				printf("%d -> received index %d\n", my_rank, received[x].index);
				a[received[x].index].real = received[x].real;
				a[received[x].index].imag = received[x].imag;
			}
		}
		*/
		cycles++;
	}
	//free(to_send);
	//free(received);
	
	/*
	if (invert) {
		int i;
		for (i = my_start; i < my_end; i++) {
			a[i].real /= n;
			a[i].imag /= n;
		}
	}*/

}

int main(int argc, char* argv[]) {
	int comm_sz;
	int my_rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Timers declaration for measuring time
	clock_t start, end;

	send_tuple data;

	// Definizione del tipo di dato MPI per send_tuple
	int blocklengths[3] = {1, 1, 1};
	MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
	MPI_Aint offsets[3];

	MPI_Aint baseaddr;
	MPI_Get_address(&data, &baseaddr);
	MPI_Get_address(&data.real, &offsets[0]);
	MPI_Get_address(&data.imag, &offsets[1]);
	MPI_Get_address(&data.index, &offsets[2]);

	offsets[0] = MPI_Aint_diff(offsets[0], baseaddr);
	offsets[1] = MPI_Aint_diff(offsets[1], baseaddr);
	offsets[2] = MPI_Aint_diff(offsets[2], baseaddr);

	MPI_Type_create_struct(3, blocklengths, offsets, types, &mpi_send_tuple_type);
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
		const char *input_file_name = "../dataset/data/dataset_0_0.txt";
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
		complex *a = malloc(n * sizeof(complex));

		// Reading input array
		int i;
		for (i = 0; i < n; i++) {
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

		for (i = 0; i < n; i++) {
			if (i < reverse(i, lg_n))
				swap(&a[i], &a[reverse(i, lg_n)]);
		}

		end = clock();

		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for reordering the input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		start = clock();

		//Broadcast n
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//Broadcast the array a
		MPI_Bcast(a, n, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		end = clock();

		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for broadcasting the input array: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		start = clock();

		//Solve my data
		parallel_fft(a, n, my_rank, comm_sz, lg_n, 0);

		end = clock();
		if (PRINTING_TIME) {
			fprintf(timings_file, "Time for calculating the parallel fft: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		}

		//Print my result


		fclose(timings_file);

		//Free memory
		free(a);
	} else {
		// Receiving n
		int n;
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Receiving the array a
		complex *a = malloc(n * sizeof(complex));
		MPI_Bcast(a, n, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		// Calculating lg_n
		int lg_n = 0;
		while ((1 << lg_n) < n)
			lg_n++;

		//Solve my data
		parallel_fft(a, n, my_rank, comm_sz, lg_n, 0);


		//Print my result


		//Free memory
		free(a);
	}

	//MPI_Type_free(&mpi_send_tuple_type);

	MPI_Finalize();
	return 0;
}
