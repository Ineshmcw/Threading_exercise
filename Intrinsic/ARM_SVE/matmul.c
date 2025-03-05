#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arm_sve.h>
#include <time.h>

#define N 128
#define M 128

// Standard matrix multiplication using 2D arrays
void mulMat(float m1[M][N], float m2[N][M], float res[M][M], int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            float sum = 0.0f;
            for (int k = 0; k < n; k++)
            {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
}


void mulMatAVX(float m1[M][N], float m2[N][M], float res[M][M], int m, int n, int vec_len)
{
    // int chunk_size = 8;
    int vl = vec_len;
    int aligned_n = (n / vl) * vl; // Largest multiple of 8 ≤ N
    int aligned_m = (m / vl) * vl; // Largest multiple of 8 ≤ N
    svbool_t pred = svptrue_b32();

    // Process full 8×8 blocks
    for (int bi = 0; bi < aligned_m; bi += vl) 
    {
        for (int bj = 0; bj < aligned_m; bj += vl)
        {
            for (int bk = 0; bk < aligned_n; bk += vl)
            {
                // Load 8 rows of m2
                for (int rowm2 = bk; rowm2 < bk + vl; rowm2 += 4)
                {
                    float debug_array[svcntw()];  // Temporary array to hold vector values

                    svfloat32_t row0 = svld1_f32(pred, &m2[rowm2 + 0][bj]);
                    svfloat32_t row1 = svld1_f32(pred, &m2[rowm2 + 1][bj]);
                    svfloat32_t row2 = svld1_f32(pred, &m2[rowm2 + 2][bj]);
                    svfloat32_t row3 = svld1_f32(pred, &m2[rowm2 + 3][bj]);

                    for (int rowm1 = bi; rowm1 < bi + vl; rowm1++)
                    {
                        svfloat32_t res_store = svld1_f32( pred, &res[rowm1][bj] );

                        res_store = svadd_f32_m(pred, res_store, svmul_f32_m(pred, svdup_n_f32(m1[rowm1][rowm2 + 0]), row0));
                        res_store = svadd_f32_m(pred, res_store, svmul_f32_m(pred, svdup_n_f32(m1[rowm1][rowm2 + 1]), row1));
                        res_store = svadd_f32_m(pred, res_store, svmul_f32_m(pred, svdup_n_f32(m1[rowm1][rowm2 + 2]), row2));
                        res_store = svadd_f32_m(pred, res_store, svmul_f32_m(pred, svdup_n_f32(m1[rowm1][rowm2 + 3]), row3));

                        svst1_f32(pred, &res[rowm1][bj], res_store);
                    }
                }
            }

            // // Load the remaining rows of m2
            // for (int rowm2 = aligned_n; rowm2 < n; rowm2++)
            // {
            //     __m256 row = _mm256_loadu_ps(&m2[rowm2][bj]);

            //     for (int rowm1 = bi; rowm1 < bi + vl; rowm1++)
            //     {
            //         __m256 res_store = _mm256_loadu_ps(&res[rowm1][bj]);
            //         res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2]), row));
            //         _mm256_storeu_ps(&res[rowm1][bj], res_store);
            //     }
            // }
        }
    }

    // // Handling remaining rows and columns
    // for (int i = 0; i < m; i++)
    // {
    //     for (int j = aligned_n; j < m; j++)
    //     { // Remaining columns
    //         float sum = 0;
    //         for (int k = 0; k < n; k++)
    //         {
    //             sum += m1[i][k] * m2[k][j];
    //         }
    //         res[i][j] = sum;
    //     }
    // }
    
    // for (int i = aligned_m; i < m; i++)
    // { // Remaining rows
    //     for (int j = 0; j < aligned_n; j++)
    //     {
    //         float sum = 0;
    //         for (int k = 0; k < n; k++)
    //         {
    //             sum += m1[i][k] * m2[k][j];
    //         }
    //         res[i][j] = sum;
    //     }
    // }
}

// void _print_m256(__m256 val)
// {
//     float *f = (float *)&val;
//     printf("{%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f}\n",
//            f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
// }

void printMatrices(float res[M][M], int isAVX)
{
    if (isAVX)
        printf("Result Matrix(after AVX):\n");
    else
        printf("Result Matrix(before AVX):\n");
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            printf("%.2f ", res[i][j]);
        }
        printf("\n");
    }
}

// void mulMatAVX(float m1[N][N], float m2[N][N], float res[N][N], int n)
// {
//     float trans[N][N] = {0};

//     // Transpose m2 into trans using SVE
//     for (int i = 0; i < N; i += svcntw()) {
//         for (int j = 0; j < N; j++) {
//             svbool_t pred = svwhilelt_b32(i, N);
//             svfloat32_t row = svld1(pred, &m2[i][j]);
//             svst1(pred, &trans[j][i], row);
//         }
//     }

//     // Matrix multiplication using SVE
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             svfloat32_t sum = svdup_f32(0.0f);

//             for (int k = 0; k < N; k += svcntw()) {
//                 svbool_t pred = svwhilelt_b32(k, N);
//                 svfloat32_t vecA = svld1(pred, &m1[i][k]);
//                 svfloat32_t vecB = svld1(pred, &trans[j][k]);
//                 sum = svmla_f32_m(pred, sum, vecA, vecB);
//             }

//             // Reduce sum across the vector
//             res[i][j] = svaddv_f32(svptrue_b32(), sum);
//         }
//     }
// }

int areMatricesEqual(float m1[M][M], float m2[M][M], float tolerance)
{
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            float diff = fabs(m1[i][j] - m2[i][j]);
            if (diff > tolerance)
            {
                printf("m1[%d][%d] = %.2f, m2[%d][%d] = %.2f\n", i, j, m1[i][j], i, j, m2[i][j]);
                printf("Difference: %.2f\n", diff);
                return 0; // Matrices are not equal
            }
        }
    }
    printf("Matrices are equal within the tolerance.\n");
    return 1; // Matrices are equal within the tolerance
}

int main()
{
    float m1[M][N] __attribute__((aligned(32)));
    float m2[N][M] __attribute__((aligned(32)));

    float res1[M][M] = {0}, res2[M][M] = {0};

    int vl = svcntw();
    printf("SVE vector length: %d\n", vl);

    printf("m1 matrix:\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            m1[i][j] = 1;
            printf("%.2f ", m1[i][j]);
        }
        printf("\n");
    }

    printf("m2 matrix:\n");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
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
    mulMatAVX(m1, m2, res2, M, N, vl);
    end = clock();
    double elapsed_avx = (double)(end - start) / CLOCKS_PER_SEC;

    printMatrices(res1, 0);
    printMatrices(res2, 1);

    printf("Execution time (Before AVX): %.6f seconds\n", elapsed_non_avx);
    printf("Execution time (With AVX): %.6f seconds\n", elapsed_avx);

    double speed_difference = ((elapsed_non_avx - elapsed_avx) / elapsed_non_avx) * 100.0;
    printf("Speed improvement with AVX: %.2f%%\n", (double)speed_difference);

    // Now check if the results are equal
    if (areMatricesEqual(res1, res2, 1e-6))
    {
        printf("The results of the two matrix multiplication functions are equal!\n");
    }
    else
    {
        printf("The results of the two matrix multiplication functions are NOT equal!\n");
    }

    return 0;
}




