#include <stdio.h>
#include <omp.h>
#include <math.h>

int is_prime(int n) {
    if (n < 2) return 0;
    for (int i = 2; i <= sqrt(n); i++) {
        if (n % i == 0) return 0;
    }
    return 1;
}

int main() {
    int n = 100000; // Check numbers from 1 to n
    int prime_count = 0;

    double start_time, run_time;

    // Start timing
    start_time = omp_get_wtime();

    #pragma omp parallel for schedule(static, 100) reduction(+: prime_count)
    for (int i = 1; i <= n; i++) {
        if (is_prime(i)) {
            prime_count++;
        }
    }

    // End timing
    run_time = omp_get_wtime() - start_time;

    printf("Number of primes between 1 and %d is %d\n", n, prime_count);
    printf("Dynamic scheduling took %f seconds\n", run_time);

    return 0;
}
