/* Include benchmark-specific header. */
#define LARGE_DATASET

#include "gemver.h"

double bench_t_start, bench_t_end;
int nProcs, rank;

void bench_timer_start() {
    bench_t_start = MPI_Wtime();
}

void bench_timer_stop() {
    bench_t_end = MPI_Wtime();
}

void bench_timer_print() {
    if (rank == 0) {
        printf("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
    }
}


int start_index, size;

static
void init_array(int n,
                float *alpha,
                float *beta,
                float A[size][n],
                float u1[n],
                float v1[n],
                float u2[n],
                float v2[n],
                float w[n],
                float x[n],
                float y[n],
                float z[n]) {
    int i, j;
    *alpha = 1.5;
    *beta = 1.2;

    float fn = (float) n;

    for (i = 0; i < n; i++) {
        u1[i] = i;
        u2[i] = ((i + 1) / fn) / 2.0;
        v1[i] = ((i + 1) / fn) / 4.0;
        v2[i] = ((i + 1) / fn) / 6.0;
        y[i] = ((i + 1) / fn) / 8.0;
        z[i] = ((i + 1) / fn) / 9.0;
        x[i] = 0.0;
        w[i] = 0.0;
    }
    for (i = 0; i < size; ++i)
        for (j = 0; j < n; j++)
            A[i][j] = (float) ((i + start_index) * j % n) / n;
}

static
void print_array(int n,
                 float w[n]) {
    if (rank != 0) {
        return;
    }
    int i;
    fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
    fprintf(stderr, "begin dump: %s", "w");
    for (i = 0; i < n; i++) {
        if (i % 20 == 0) fprintf(stderr, "\n");
        fprintf(stderr, "%0.2f ", w[i]);
    }
    fprintf(stderr, "\nend   dump: %s\n", "w");
    fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

static
void kernel_gemver(int n,
                   float alpha,
                   float beta,
                   float A[size][n],
                   float u1[n],
                   float v1[n],
                   float u2[n],
                   float v2[n],
                   float w[n],
                   float x[n],
                   float y[n],
                   float z[n]) {
    int i, j;

    for (i = 0; i < size; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = A[i][j] + u1[i + start_index] * v1[j] + u2[i + start_index] * v2[j];
        }
    }

    for (i = 0; i < size; i++) {
        for (j = 0; j < n; j++) {
            x[j] = x[j] + beta * A[i][j] * y[i + start_index];
        }
    }

    MPI_Status status;
    if (rank == 0) {
        float (*temp_x)[n];
        temp_x = (float (*)[n]) malloc((n) * sizeof(float));

        for (i = 0; i < n; i++) {
            x[i] = x[i] + z[i];
        }

        for (i = 1; i < nProcs; ++i) {
            MPI_Recv(*temp_x, n, MPI_FLOAT, i, 23, MPI_COMM_WORLD, &status);
            for (j = 0; j < n; j++) {
                x[j] = x[j] + (*temp_x)[j];
            }
        }

        for (i = 1; i < nProcs; i++) {
            MPI_Send(x, n, MPI_FLOAT, i, 2312, MPI_COMM_WORLD);
        }

        free((void *) temp_x);

    } else {
        MPI_Send(x, n, MPI_FLOAT, 0, 23, MPI_COMM_WORLD);
        MPI_Recv(x, n, MPI_FLOAT, 0, 2312, MPI_COMM_WORLD, &status);

    }

    for (i = 0; i < size; i++) {
        for (j = 0; j < n; j++) {
            w[i + start_index] = w[i + start_index] + alpha * A[i][j] * x[j];
        }
    }

}

int min(int a, int b) {
    return (a < b) ? a : b;
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = N;
    start_index = n / nProcs * rank + min(rank, n % nProcs);
    size = n / nProcs + (rank < n % nProcs);

    float alpha;
    float beta;
    float (*A)[size][n];
    A = (float (*)[size][n]) malloc((size) * (n) * sizeof(float));
    float (*u1)[n];
    u1 = (float (*)[n]) malloc((n) * sizeof(float));
    float (*v1)[n];
    v1 = (float (*)[n]) malloc((n) * sizeof(float));
    float (*u2)[n];
    u2 = (float (*)[n]) malloc((n) * sizeof(float));
    float (*v2)[n];
    v2 = (float (*)[n]) malloc((n) * sizeof(float));
    float (*w)[n];
    w = (float (*)[n]) malloc((n) * sizeof(float));
    float (*x)[n];
    x = (float (*)[n]) malloc((n) * sizeof(float));
    float (*y)[n];
    y = (float (*)[n]) malloc((n) * sizeof(float));
    float (*z)[n];
    z = (float (*)[n]) malloc((n) * sizeof(float));

    init_array(n, &alpha, &beta,
               *A,
               *u1,
               *v1,
               *u2,
               *v2,
               *w,
               *x,
               *y,
               *z);

    bench_timer_start();

    kernel_gemver(n, alpha, beta,
                  *A,
                  *u1,
                  *v1,
                  *u2,
                  *v2,
                  *w,
                  *x,
                  *y,
                  *z);

    if (rank == 0) {
        MPI_Status status;
        int j;
        for (j = 1; j < nProcs; j++) {
            int temp_size = n / nProcs + (j < n % nProcs);
            MPI_Recv(&((*w)[n / nProcs * j + min(j, n % nProcs)]), temp_size, MPI_FLOAT, j, 231299, MPI_COMM_WORLD,
                     &status);
        }
    } else {
        MPI_Send(&((*w)[start_index]), size, MPI_FLOAT, 0, 231299, MPI_COMM_WORLD);
    }

    bench_timer_stop();

    bench_timer_print();

//    print_array(n, *w);

    free((void *) A);
    free((void *) u1);
    free((void *) v1);
    free((void *) u2);
    free((void *) v2);
    free((void *) w);
    free((void *) x);
    free((void *) y);
    free((void *) z);

    MPI_Finalize();

    return 0;
}
