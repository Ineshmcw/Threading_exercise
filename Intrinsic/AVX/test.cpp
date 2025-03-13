#include <iostream>
#include <arm_sve.h>
#include <chrono>

#define N 512  // Increase matrix size for meaningful performance comparison
#define ITERATIONS 1000  // Number of runs for timing

void transpose_sve(float m2[N][N], float trans[N][N]) {
    for (int rows = 0; rows < N; rows++) {
        for (int columns = 0; columns < N; columns += svcntw()) {
            svbool_t pg = svwhilelt_b32(0, N - columns);
            svfloat32_t row = svld1(pg, &m2[rows][columns]);

            // Store elements manually
            float temp[svcntw()];
            svst1(pg, temp, row);

            for (int k = 0; k < svcntw() && (columns + k) < N; k++) {
                trans[columns + k][rows] = temp[k];
            }
        }
    }
}

void transpose_normal(float m2[N][N], float trans[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            trans[j][i] = m2[i][j];
        }
    }
}

void benchmark() {
    float m2[N][N], trans[N][N];

    // Initialize the matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            m2[i][j] = i * N + j + 1;
        }
    }

    // Measure SVE transpose
    auto start_sve = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ITERATIONS; i++) {
        transpose_sve(m2, trans);
    }
    auto end_sve = std::chrono::high_resolution_clock::now();
    double time_sve = std::chrono::duration<double, std::milli>(end_sve - start_sve).count();

    // Measure Normal transpose
    auto start_normal = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ITERATIONS; i++) {
        transpose_normal(m2, trans);
    }
    auto end_normal = std::chrono::high_resolution_clock::now();
    double time_normal = std::chrono::duration<double, std::milli>(end_normal - start_normal).count();

    // Print timing results
    std::cout << "SVE Transpose Time: " << time_sve << " ms\n";
    std::cout << "Normal Transpose Time: " << time_normal << " ms\n";
    std::cout << "Speedup Factor: " << time_normal / time_sve << "x\n";
}

int main() {
    benchmark();
    return 0;
}
