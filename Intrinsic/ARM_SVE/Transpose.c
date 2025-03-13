/*
 * To run:
aarch64-linux-gnu-gcc -g -O2 -march=armv8-a+sve -o Transpose Transpose.c
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./transpose
 */

#include <stdio.h>
#include <arm_sve.h>


#define ROWS 8
#define COLS 8

void transpose(float A[ROWS][COLS], float result[ROWS][COLS]) {
    svbool_t pred = svptrue_b32();
    // Transpose m2 into trans using SVE
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j += svcntw()) {
            svfloat32_t row = svld1_f32(pred, &A[i][j]);
            svst1(pred, &result[j][i], row);


            printf("%5.1f",A[i][j]);
            printf("%5.1f\n",result[j][i]);
        }
        printf("\n");
    }
}
void printMatrix(float matrix[ROWS][COLS]) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            printf("%5.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    float A[ROWS][COLS] = {0};

    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            A[i][j] = i+j*(2+i);
        }
    }


    float result[ROWS][COLS] = {0};

    transpose(A, result);

    printf("Matrix A:\n");
    printMatrix(A);

    printf("Transpose of Matrix:\n");
    printMatrix(result);

    return 0;
}

