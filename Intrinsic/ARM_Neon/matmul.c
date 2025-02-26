#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <arm_neon.h>

#define N 128
#define M 128

// To run
// arm-linux-gnueabihf-gcc -march=armv7-a -mfpu=neon -mfloat-abi=hard -o matmul matmul.c
// qemu-arm -L /usr/arm-linux-gnueabihf/ ./matmul

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


void mulMatNEON(float m1[M][N], float m2[N][M], float res[M][M], int m, int n) {
    float trans[N][N] __attribute__((aligned(32))) = {0};

    for (int i = 0; i < N; i += 4) {
        for (int j = 0; j < N; j += 4) {
            float32x4_t row0 = vld1q_f32(&m2[i + 0][j]);
            float32x4_t row1 = vld1q_f32(&m2[i + 1][j]);
            float32x4_t row2 = vld1q_f32(&m2[i + 2][j]);
            float32x4_t row3 = vld1q_f32(&m2[i + 3][j]);

            // **Transpose 4x4 using NEON**
            float32x4x2_t t01 = vzipq_f32(row0, row1);
            float32x4x2_t t23 = vzipq_f32(row2, row3);

            float32x4_t tt0 = vcombine_f32(vget_low_f32(t01.val[0]), vget_low_f32(t23.val[0]));
            float32x4_t tt1 = vcombine_f32(vget_high_f32(t01.val[0]), vget_high_f32(t23.val[0]));
            float32x4_t tt2 = vcombine_f32(vget_low_f32(t01.val[1]), vget_low_f32(t23.val[1]));
            float32x4_t tt3 = vcombine_f32(vget_high_f32(t01.val[1]), vget_high_f32(t23.val[1]));

            vst1q_f32(&trans[j + 0][i], tt0);
            vst1q_f32(&trans[j + 1][i], tt1);
            vst1q_f32(&trans[j + 2][i], tt2);
            vst1q_f32(&trans[j + 3][i], tt3);
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            float32x4_t sum = vdupq_n_f32(0);

            for (int k = 0; k < N; k += 4) {
                float32x4_t vecA = vld1q_f32(&m1[i][k]);
                float32x4_t vecB = vld1q_f32(&trans[j][k]);

                sum = vfmaq_f32(sum, vecA, vecB);
            }

            // **Horizontal Sum**
            float32x2_t sum2 = vadd_f32(vget_low_f32(sum), vget_high_f32(sum));
            float final_sum = vget_lane_f32(vpadd_f32(sum2, sum2), 0);

            res[i][j] = final_sum;
        }
    }
}

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
    mulMatNEON(m1, m2, res2, M, N);
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
