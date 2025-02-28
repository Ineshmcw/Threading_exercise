#include <stdio.h>
#include <arm_sve.h>


#define ROWS 128
#define COLS 128

void addMatrices(int A[ROWS][COLS], int B[ROWS][COLS], int result[ROWS][COLS]) {
    svbool_t pred;
    for (int i = 0; i < COLS; i+= svcntw()) {
        svbool_t pred = svwhilelt_b32(i, COLS);
        for (int j = 0; j < ROWS; j++) {
            svint32_t vec_a = svld1_s32(pred, &A[j][i]);
            svint32_t vec_b = svld1_s32(pred, &B[j][i]);

            svint32_t vec_sum = svadd_s32_m(pred, vec_a, vec_b);
            svst1_s32(pred, &result[j][i], vec_sum);
        }
    }
}

void printMatrix(int matrix[ROWS][COLS]) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int A[ROWS][COLS] = {0};
    int B[ROWS][COLS] = {0};

    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            A[i][j] = 1;
            // printf("%.2f ", m1[i][j]);
        }
        // printf("\n");
    }

    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            B[i][j] = 2;
            // printf("%.2f ", m2[i][j]);
        }
        // printf("\n");
    }


    int result[ROWS][COLS];

    addMatrices(A, B, result);

    printf("Matrix A:\n");
    printMatrix(A);

    printf("\nMatrix B:\n");
    printMatrix(B);

    printf("\nSum of matrices:\n");
    printMatrix(result);

    return 0;
}