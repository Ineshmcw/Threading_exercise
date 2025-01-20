#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

static long num_steps = 100000000;
double step;

int main()
{
    int i;
    double pi = 0.0;
    double start_time, run_time;

    step = 1.0 / (double)num_steps;
    int no_threads = 4;
    double *thread_sum = malloc(no_threads * sizeof(double));

    if (thread_sum == NULL)
    {
        printf("Memory allocation failed\n");
        return 1;
    }

    start_time = omp_get_wtime();

    #pragma omp parallel num_threads(no_threads)
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        int start = (thread_id * num_steps) / num_threads;
        int end = ((thread_id + 1) * num_steps) / num_threads;

        double local_sum = 0.0;

        for (int i = start; i < end; i++)
        {
            double x = (i + 0.5) * step;
            local_sum += 4.0 / (1.0 + x * x);
        }
        thread_sum[thread_id] = local_sum;
        printf("Thread %d sum %lf\n", thread_id, local_sum);
    }

    for (i = 0; i < no_threads; i++)
    {
        pi += thread_sum[i];
    }
    pi *= step;

    run_time = omp_get_wtime() - start_time;
    printf("\n pi with %ld steps is %lf in %lf seconds\n", num_steps, pi, run_time);

    free(thread_sum); // Free allocated memory

    return 0;
}
