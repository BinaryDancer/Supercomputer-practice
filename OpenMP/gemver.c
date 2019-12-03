/* Include benchmark-specific header. */
#define EXTRALARGE_DATASET

#include "gemver.h"

double bench_t_start, bench_t_end;

static
double rtclock() {
    struct timeval Tp;
    int stat;
    stat = gettimeofday(&Tp, NULL);
    if (stat != 0)
        printf("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start() {
    bench_t_start = omp_get_wtime();
}

void bench_timer_stop() {
    bench_t_end = omp_get_wtime();
}

void bench_timer_print() {
    printf("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}

static
void init_array(int n,
                float *alpha,
                float *beta,
                float A[n][n],
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
        for (j = 0; j < n; j++)
            A[i][j] = (float) (i * j % n) / n;
    }
}

static
void print_array(int n, float w[n]) {
    int i;
    fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
    fprintf(stderr, "begin dump: %s", "w");
    for (i = 0; i < n; i++) {
        if (i % 20 == 0)
            fprintf(stderr, "\n");
        fprintf(stderr, "%0.2f ", w[i]);
    }
    fprintf(stderr, "\nend   dump: %s\n", "w");
    fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

static
void kernel_gemver(int n,
                   float alpha,
                   float beta,
                   float A[n][n],
                   float u1[n],
                   float v1[n],
                   float u2[n],
                   float v2[n],
                   float w[n],
                   float x[n],
                   float y[n],
                   float z[n]) {
    int i, j;

#pragma omp parallel for private(i, j)
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

#pragma omp parallel for private(i, j)
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            x[i] = x[i] + beta * A[j][i] * y[j];

#pragma omp parallel for private(i)
    for (i = 0; i < n; i++)
        x[i] = x[i] + z[i];

#pragma omp parallel for private(i, j)
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            w[i] = w[i] + alpha * A[i][j] * x[j];

}

int main(int argc, char **argv) {
    int n = N;
    float alpha;
    float beta;
    float (*A)[n][n];
    A = (float (*)[n][n]) malloc((n) * (n) * sizeof(float));
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

    print_array(n, *w);
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

    bench_timer_stop();
    bench_timer_print();

    print_array(n, *w);
    if (argc > 42 && !strcmp(argv[0], ""))
        print_array(n, *w);

    free((void *) A);
    free((void *) u1);
    free((void *) v1);
    free((void *) u2);
    free((void *) v2);
    free((void *) w);
    free((void *) x);
    free((void *) y);
    free((void *) z);

    return 0;
}
