#include <stdio.h>
#include <omp.h>
static long num_steps = 100000000;
double step;

int main()
{
    int i;
    double x, pi = 0.0;
    double start_time, run_time;

    step = 1.0 / (double)num_steps;
    start_time = omp_get_wtime();

    double thread_sum[2] = {0.0, 0.0};

#pragma omp parallel num_threads(2)
    {
        int thread_id = omp_get_thread_num();
        int start = (thread_id * num_steps) / 2 + 1;
        int end = ((thread_id + 1) * num_steps) / 2;

        for (int i = start; i <= end; i++)
        {
            double x = (i - 0.5) * step;
            thread_sum[thread_id] += 4.0 / (1.0 + x * x);
        }
        printf("Thread %d sum %lf\n", thread_id, thread_sum[thread_id]);
    }

    for (i = 0; i < 2; i++)
    {
        pi += thread_sum[i];
    }
    pi *= step;

    run_time = omp_get_wtime() - start_time;
    printf("\n pi with %ld steps is %lf in %lf seconds\n", num_steps, pi, run_time);

    return 0;
}
