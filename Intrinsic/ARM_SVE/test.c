#include <arm_sve.h>
#include <stdio.h>

#define N 8
#define M 8

void mulMatAVX(int vec_len)
{
    int vl = svcntw();

    printf("function in%d\n",vec_len);
    printf("vector len inside fun:%d",vl);
}

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


int main() {
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
    int vl = svcntw();
    printf("SVE vector length: %d\n", vl);
    mulMatAVX(vl);

    printMatrices(res2, 1);


    return 0;
}
