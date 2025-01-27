#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// Function to calculate (base^exp) % mod
uint64_t modular_exponentiation(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base = base % mod; 
    if (base == 0) return 0;

    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}

// Function to calculate modular exponentiation in chunks
uint64_t parallel_modular_exponentiation(uint64_t base, uint64_t exp, uint64_t mod) {
    const uint64_t CHUNK_SIZE = 1000000000;
    uint64_t partial_result = 1;

    while (exp > 0) {
        uint64_t current_chunk = (exp > CHUNK_SIZE) ? CHUNK_SIZE : exp;
        uint64_t chunk_result = modular_exponentiation(base, current_chunk, mod);
        partial_result = (partial_result * chunk_result) % mod;
        exp -= current_chunk;
    }

    return partial_result;
}

int main() {
    srand(time(0));
    int num_tests = 5;  // Number of test cases

    printf("Generating %d test cases with larger inputs...\n\n", num_tests);

    for (int i = 1; i <= num_tests; i++) {
        // Generate larger random base, exponent, and modulus
        uint64_t base = (uint64_t)rand() * rand() % 1000000000 + 1;  // Base up to 10^9
        uint64_t exp = (uint64_t)rand() * rand() % 1000000000 + 1;   // Exponent up to 10^9
        uint64_t mod = (uint64_t)rand() * rand() % 1000000000 + 2;   // Modulus up to 10^9

        double start_time = omp_get_wtime();
        uint64_t result = parallel_modular_exponentiation(base, exp, mod);
        double end_time = omp_get_wtime();
        printf("  Time Taken: %f seconds\n\n", end_time - start_time);
    }

    return 0;
}
