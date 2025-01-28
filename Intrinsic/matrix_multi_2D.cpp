#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <math.h>

using namespace std;

const int N = 8;

// Standard matrix multiplication using 2D arrays
void mulMat(float m1[N][N], float m2[N][N], float res[N][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float sum = 0.0f;
            for (int k = 0; k < n; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
}

// void mulMatAVX(float m1[N][N], float m2[N][N], float res[N][N], int n) {
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             __m256 sum = _mm256_setzero_ps();
//             for (int k = 0; k < n; k += 8) {
//                 __m256 row = _mm256_loadu_ps(&m1[i][k]);
//                 __m256 col = _mm256_loadu_ps(&m2[k][j]);

//                 sum = _mm256_add_ps(sum, _mm256_mul_ps(row, col));
//             }

//             float temp[8];
//             _mm256_storeu_ps(temp, sum);
//             float sum_s = 0.0f;
//             for (int k = 0; k < 8; ++k)
//                 sum_s += temp[k];

//             res[i][j] = sum_s;
//         }
//     }
// }

void mulMatAVX(float m1[N][N], float m2[N][N], float res[N][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            __m256 sum = _mm256_setzero_ps();

            for (int k = 0; k < n; k += 8) {
                __m256 row = _mm256_loadu_ps(&m1[i][k]);
                
                alignas(32) float colArray[8];
                for (int idx = 0; idx < 8; ++idx) {
                    colArray[idx] = m2[k + idx][j];
                }
                __m256 col = _mm256_loadu_ps(colArray);

                sum = _mm256_add_ps(sum, _mm256_mul_ps(row, col));
            }

            float temp[8];
            _mm256_storeu_ps(temp, sum);
            float sum_s = 0.0f;
            for (int k = 0; k < 8; ++k)
                sum_s += temp[k];

            res[i][j] = sum_s;
        }
    }
}

// Function to check if two matrices are equal element-wise
bool areMatricesEqual(float m1[N][N], float m2[N][N], float tolerance = 1e-6) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Compare the absolute difference with tolerance
            if (fabs(m1[i][j] - m2[i][j]) > tolerance) {
                return false; // Matrices are not equal
            }
            // else {
            //     cout << "m1[" << i << "][" << j << "]: " << m1[i][j] << " m2[" << i << "][" << j << "]: " << m2[i][j] << endl;
            // }
        }
    }
    return true; // Matrices are equal within the tolerance
}

int main() {
    float m1[N][N], m2[N][N], res1[N][N] = {0}, res2[N][N] = {0};

    // Initialize matrices with 1s
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            m1[i][j] = 1.0f; // Fill m1 with 1s
            m2[i][j] = 1.0f; // Fill m2 with 1s
        }
    }

    // Measure time for non-AVX implementation using 2D arrays
    auto start = chrono::high_resolution_clock::now();
    mulMat(m1, m2, res1, N);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_non_avx = end - start;

    cout << "Execution time (Before AVX): " << elapsed_non_avx.count() << " seconds\n";

    // Measure time for AVX-based implementation
    start = chrono::high_resolution_clock::now();
    mulMatAVX(m1, m2, res2, N);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_avx = end - start;

    cout << "Execution time (With AVX): " << elapsed_avx.count() << " seconds\n";

    // Now check if the results are equal
    if (areMatricesEqual(res1, res2)) {
        cout << "The results of the two matrix multiplication functions are equal!" << endl;
    } else {
        cout << "The results of the two matrix multiplication functions are NOT equal!" << endl;
    }

    return 0;
}