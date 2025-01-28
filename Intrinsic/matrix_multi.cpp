#include <iostream>
#include <chrono>
#include <immintrin.h>

using namespace std;

const int N = 512; // Size of the square matrix

void mulMat(float* m1, float* m2, float* res, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float sum = 0.0f;
            for (int k = 0; k < n; k++) {
                sum += m1[i * n + k] * m2[k * n + j];
            }
            res[i * n + j] = sum;
        }
    }
}

void mulMatAVX(float* m1, float* m2, float* res, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            __m256 sum = _mm256_setzero_ps();
            for (int k = 0; k < n; k += 8) {
                __m256 row = _mm256_loadu_ps(&m1[i * n + k]);
                __m256 col = _mm256_loadu_ps(&m2[k * n + j]);
                sum = _mm256_fmadd_ps(row, col, sum);
            }

            float temp[8];
            _mm256_storeu_ps(temp, sum);
            res[i * n + j] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
        }
    }
}

int main() {
    alignas(32) float m1[N * N], m2[N * N], res[N * N] = {0}; // Matrices are aligned for AVX

    // Initialize matrices
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            m1[i * N + j] = 1.0f; // Fill with 1s
            m2[i * N + j] = 1.0f; // Fill with 1s
        }
    }

    // Measure time for non-AVX implementation
    auto start = chrono::high_resolution_clock::now();
    mulMat(m1, m2, res, N);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_non_avx = end - start;

    cout << "Execution time (Before AVX): " << elapsed_non_avx.count() << " seconds\n";

    // Measure time for AVX implementation
    start = chrono::high_resolution_clock::now();
    mulMatAVX(m1, m2, res, N);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_avx = end - start;

    cout << "Execution time (After AVX): " << elapsed_avx.count() << " seconds\n";

    // Calculate speed difference in percentage
    double speed_difference = ((elapsed_non_avx.count() - elapsed_avx.count()) / elapsed_non_avx.count()) * 100.0;
    cout << "Speed improvement with AVX: " << speed_difference << "%\n";

    return 0;
}
