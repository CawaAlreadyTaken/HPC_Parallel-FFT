#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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

void fft(complex *a, int n, int invert) {
    int lg_n = 0;
    while ((1 << lg_n) < n)
        lg_n++;

    for (int i = 0; i < n; i++) {
        if (i < reverse(i, lg_n))
            swap(&a[i], &a[reverse(i, lg_n)]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        complex wlen = complex_from_polar(1.0, ang);
        for (int i = 0; i < n; i += len) {
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

    if (invert) {
        for (int i = 0; i < n; i++) {
            a[i].real /= n;
            a[i].imag /= n;
        }
    }
}

int main() {
    // Timers declaration for measuring time
    clock_t start, end;

    // Opening file for writing time results
    FILE *file = fopen("timing_serial_solver_0.txt", "w");

    // Reading input size
    int n;
    scanf("%d", &n);

    start = clock();

    // Allocating memory for input array
    complex *a = malloc(n * sizeof(complex));

    // Reading input array
    for (int i = 0; i < n; i++) {
        scanf("%lf", &a[i].real);
        a[i].imag = 0;
    }

    end = clock();
    if (PRINTING_TIME) {
        fprintf(file, "Time for allocating memory and reading input: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    }

    start = clock();

    // Calling serial FFT function
    fft(a, n, 0);

    end = clock();
    if (PRINTING_TIME) {
        fprintf(file, "Time for FFT: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    }

    if (PRINTING_OUTPUT) {
        start = clock();

        printf("FFT result:\n");
        for (int i = 0; i < n; i++) {
            printf("(%f, %f)\n", a[i].real, a[i].imag);
        }

        end = clock();
        if (PRINTING_TIME) {
            fprintf(file, "Time for printing output: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        }
    }

    // Close file
    fclose(file);

    // Freeing memory
    free(a);
    return 0;
}

