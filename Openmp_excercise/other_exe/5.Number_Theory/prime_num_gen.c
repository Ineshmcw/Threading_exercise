// C program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

void SieveOfEratosthenes(int n)
{
  
    // Create a boolean array "prime[0..n]" and initialize
    // all entries it as true. A value in prime[i] will
    // finally be false if i is Not a prime, else true.
    bool prime[n + 1];
    memset(prime, true, sizeof(prime));

    for (int p = 2; p * p <= n; p++) {
        // If prime[p] is not changed, then it is a prime
        if (prime[p] == true) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked.
            for (int i = p * p; i <= n; i += p)
                prime[i] = false;
        }
    }

}

// Driver Code
int main() {
    int n;
    n = 1000000;

    if (n < 2) {
        return 0;
    }

    // Measure execution time
    double start_time = omp_get_wtime();
    SieveOfEratosthenes(n);
    double end_time = omp_get_wtime();

    printf("Execution time: %f seconds\n", end_time - start_time);
    return 0;
}