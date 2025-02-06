#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>
#include <time.h>

#define N 256
#define M 128

// Standard matrix multiplication using 2D arrays
void mulMat(float m1[M][N], float m2[N][M], float res[M][M], int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            float sum = 0.0f;
            for (int k = 0; k < n; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
}

void _print_m256(__m256 val) {
    float* f = (float*)&val;
    printf("{%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f}\n", 
           f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
}

void printMatrices(float res[M][M], int isAVX) {
    if (isAVX)
        printf("Result Matrix(after AVX):\n");
    else
        printf("Result Matrix(before AVX):\n");
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            printf("%.2f ", res[i][j]);
        }
        printf("\n");
    }
}

void mulMatAVX(float m1[M][N], float m2[N][M], float res[M][M], int m, int n) {
    int chunk_size = 8;
    int aligned_n = (n / chunk_size) * chunk_size;  // Largest multiple of 8 ≤ N
    int aligned_m = (m / chunk_size) * chunk_size;  // Largest multiple of 8 ≤ N

    // Process full 8×8 blocks
    for (int bi = 0; bi < aligned_m; bi += chunk_size) {
        for (int bj = 0; bj < aligned_m; bj += chunk_size) {
            for (int bk = 0; bk < aligned_n; bk += chunk_size) {
                // Load 8 rows of m2
                for (int rowm2 = bk; rowm2 < bk + chunk_size; rowm2 += 8) {
                    __m256 row0 = _mm256_loadu_ps(&m2[rowm2 + 0][bj]);
                    __m256 row1 = _mm256_loadu_ps(&m2[rowm2 + 1][bj]);
                    __m256 row2 = _mm256_loadu_ps(&m2[rowm2 + 2][bj]);
                    __m256 row3 = _mm256_loadu_ps(&m2[rowm2 + 3][bj]);
                    __m256 row4 = _mm256_loadu_ps(&m2[rowm2 + 4][bj]);
                    __m256 row5 = _mm256_loadu_ps(&m2[rowm2 + 5][bj]);
                    __m256 row6 = _mm256_loadu_ps(&m2[rowm2 + 6][bj]);
                    __m256 row7 = _mm256_loadu_ps(&m2[rowm2 + 7][bj]);

                    for (int rowm1 = bi; rowm1 < bi + chunk_size; rowm1++) {
                        __m256 res_store = _mm256_loadu_ps(&res[rowm1][bj]);

                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 0]), row0));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 1]), row1));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 2]), row2));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 3]), row3));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 4]), row4));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 5]), row5));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 6]), row6));
                        res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2 + 7]), row7));

                        _mm256_storeu_ps(&res[rowm1][bj], res_store);
                    }
                }
            }

            // Load the remaining rows of m2
            for (int rowm2 = aligned_n; rowm2 < n; rowm2++) {
                __m256 row = _mm256_loadu_ps(&m2[rowm2][bj]);

                for (int rowm1 = bi; rowm1 < bi + chunk_size; rowm1++) {
                    __m256 res_store = _mm256_loadu_ps(&res[rowm1][bj]);
                    res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2]), row));
                    _mm256_storeu_ps(&res[rowm1][bj], res_store);
                }
            }
        }
    }

    // Handling remaining rows and columns
    for (int i = 0; i < n; i++) {
        for (int j = aligned_n; j < n; j++) {  // Remaining columns
            float sum = 0;
            for (int k = 0; k < n; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
    for (int i = aligned_n; i < n; i++) {  // Remaining rows
        for (int j = 0; j < aligned_n; j++) {
            float sum = 0;
            for (int k = 0; k < n; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
}

int areMatricesEqual(float m1[M][M], float m2[M][M], float tolerance) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            float diff = fabs(m1[i][j] - m2[i][j]);
            if (diff > tolerance) {
                printf("m1[%d][%d] = %.2f, m2[%d][%d] = %.2f\n", i, j, m1[i][j], i, j, m2[i][j]);
                printf("Difference: %.2f\n", diff);
                return 0;  // Matrices are not equal
            }
        }
    }
    printf("Matrices are equal within the tolerance.\n");
    return 1;  // Matrices are equal within the tolerance
}

int main() {
    float m1[M][N] __attribute__((aligned(32)));
    float m2[N][M] __attribute__((aligned(32)));

    float res1[M][M] = {0}, res2[M][M] = {0};

    printf("m1 matrix:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            m1[i][j] = 1;
            printf("%.2f ", m1[i][j]);
        }
        printf("\n");
    }

    printf("m2 matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            m2[i][j] = 2;
            printf("%.2f ", m2[i][j]);
        }
        printf("\n");
    }

    // Measure time for non-AVX implementation using 2D arrays
    clock_t start = clock();
    mulMat(m1, m2, res1, M, N);
    clock_t end = clock();
    double elapsed_non_avx = (double)(end - start) / CLOCKS_PER_SEC;

    // Measure time for AVX-based implementation
    start = clock();
    mulMatAVX(m1, m2, res2, M, N);
    end = clock();
    double elapsed_avx = (double)(end - start) / CLOCKS_PER_SEC;

    printMatrices(res1, 0);
    printMatrices(res2, 1);

    printf("Execution time (Before AVX): %.6f seconds\n", elapsed_non_avx);
    printf("Execution time (With AVX): %.6f seconds\n", elapsed_avx);

    double speed_difference = ((elapsed_non_avx - elapsed_avx) / elapsed_non_avx) * 100.0;
    printf("Speed improvement with AVX: %.2f%%\n", (double)speed_difference);


    // Now check if the results are equal
    if (areMatricesEqual(res1, res2, 1e-6)) {
        printf("The results of the two matrix multiplication functions are equal!\n");
    } else {
        printf("The results of the two matrix multiplication functions are NOT equal!\n");
    }

    return 0;
}
