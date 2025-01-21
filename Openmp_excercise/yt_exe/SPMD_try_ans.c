#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000;
double step;
#define PAD 8
#define NUM_THREADS 4

int main()
{
    int i, nthreads;
    double pi, sum[NUM_THREADS][PAD];
    double start_time, run_time;

    step = 1.0 / (double)num_steps;
    omp_set_num_threads(NUM_THREADS);
    start_time = omp_get_wtime();

    #pragma omp parallel
    {
        int id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        if (id == 0) nthreads = nthrds; // Set nthreads based on parallel region

        for (int i = id; i < num_steps; i += nthrds) // Fix loop increment
        {
            x = (i + 0.5) * step; // Fix calculation of x
            sum[id][0] += 4.0 / (1.0 + x * x);
        }
    }

    for (i = 0, pi = 0.0; i < nthreads; i++) // Use nthreads instead of NUM_THREADS
    {
        pi += sum[i][0] * step;
    }

    run_time = omp_get_wtime() - start_time;
    printf("\n pi with %ld steps is %lf in %lf seconds\n", num_steps, pi, run_time);

    return 0;
}
