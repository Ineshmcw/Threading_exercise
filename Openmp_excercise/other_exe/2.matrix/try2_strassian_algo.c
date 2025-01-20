/**
 * @file
 * @brief Implementation of [Strassen's
 * algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm) for matrix multiplication
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/**
 * @addtogroup matrix_algorithms Matrix algorithms
 * @{
 */

/** Allocate memory for a matrix of size n x n */
int **allocate_matrix(int n)
{
    int **matrix = (int **)malloc(n * sizeof(int *));
    if (matrix == NULL)
    {
        printf("Can't allocate memory! Please try again.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (int *)malloc(n * sizeof(int));
        if (matrix[i] == NULL)
        {
            printf("Can't allocate memory! Please try again.\n");
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

/** Free memory allocated for a matrix */
void free_matrix(int **matrix, int n)
{
    for (int i = 0; i < n; i++)
        free(matrix[i]);
    free(matrix);
}

/** Add two matrices */
void add_matrices(int **a, int **b, int **result, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i][j] = a[i][j] + b[i][j];
}

/** Subtract two matrices */
void subtract_matrices(int **a, int **b, int **result, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i][j] = a[i][j] - b[i][j];
}

/** Strassen's algorithm for matrix multiplication */
void strassen(int **a, int **b, int **result, int n)
{
    if (n == 1)
    {
        result[0][0] = a[0][0] * b[0][0];
        return;
    }

    int new_size = n / 2;
    int **a11 = allocate_matrix(new_size);
    int **a12 = allocate_matrix(new_size);
    int **a21 = allocate_matrix(new_size);
    int **a22 = allocate_matrix(new_size);
    int **b11 = allocate_matrix(new_size);
    int **b12 = allocate_matrix(new_size);
    int **b21 = allocate_matrix(new_size);
    int **b22 = allocate_matrix(new_size);
    int **c11 = allocate_matrix(new_size);
    int **c12 = allocate_matrix(new_size);
    int **c21 = allocate_matrix(new_size);
    int **c22 = allocate_matrix(new_size);
    int **m1 = allocate_matrix(new_size);
    int **m2 = allocate_matrix(new_size);
    int **m3 = allocate_matrix(new_size);
    int **m4 = allocate_matrix(new_size);
    int **m5 = allocate_matrix(new_size);
    int **m6 = allocate_matrix(new_size);
    int **m7 = allocate_matrix(new_size);
    int **temp1 = allocate_matrix(new_size);
    int **temp2 = allocate_matrix(new_size);

    // Divide matrices into quadrants
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            a11[i][j] = a[i][j];
            a12[i][j] = a[i][j + new_size];
            a21[i][j] = a[i + new_size][j];
            a22[i][j] = a[i + new_size][j + new_size];

            b11[i][j] = b[i][j];
            b12[i][j] = b[i][j + new_size];
            b21[i][j] = b[i + new_size][j];
            b22[i][j] = b[i + new_size][j + new_size];
        }
    }

    // Calculate M1 to M7
    add_matrices(a11, a22, temp1, new_size);
    add_matrices(b11, b22, temp2, new_size);
    strassen(temp1, temp2, m1, new_size);

    add_matrices(a21, a22, temp1, new_size);
    strassen(temp1, b11, m2, new_size);

    subtract_matrices(b12, b22, temp1, new_size);
    strassen(a11, temp1, m3, new_size);

    subtract_matrices(b21, b11, temp1, new_size);
    strassen(a22, temp1, m4, new_size);

    add_matrices(a11, a12, temp1, new_size);
    strassen(temp1, b22, m5, new_size);

    subtract_matrices(a21, a11, temp1, new_size);
    add_matrices(b11, b12, temp2, new_size);
    strassen(temp1, temp2, m6, new_size);

    subtract_matrices(a12, a22, temp1, new_size);
    add_matrices(b21, b22, temp2, new_size);
    strassen(temp1, temp2, m7, new_size);

    // Combine results into C
    add_matrices(m1, m4, temp1, new_size);
    subtract_matrices(temp1, m5, temp2, new_size);
    add_matrices(temp2, m7, c11, new_size);

    add_matrices(m3, m5, c12, new_size);
    add_matrices(m2, m4, c21, new_size);

    add_matrices(m1, m3, temp1, new_size);
    subtract_matrices(temp1, m2, temp2, new_size);
    add_matrices(temp2, m6, c22, new_size);

    // Combine quadrants into result
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            result[i][j] = c11[i][j];
            result[i][j + new_size] = c12[i][j];
            result[i + new_size][j] = c21[i][j];
            result[i + new_size][j + new_size] = c22[i][j];
        }
    }

    // Free allocated memory
    free_matrix(a11, new_size);
    free_matrix(a12, new_size);
    free_matrix(a21, new_size);
    free_matrix(a22, new_size);
    free_matrix(b11, new_size);
    free_matrix(b12, new_size);
    free_matrix(b21, new_size);
    free_matrix(b22, new_size);
    free_matrix(c11, new_size);
    free_matrix(c12, new_size);
    free_matrix(c21, new_size);
    free_matrix(c22, new_size);
    free_matrix(m1, new_size);
    free_matrix(m2, new_size);
    free_matrix(m3, new_size);
    free_matrix(m4, new_size);
    free_matrix(m5, new_size);
    free_matrix(m6, new_size);
    free_matrix(m7, new_size);
    free_matrix(temp1, new_size);
    free_matrix(temp2, new_size);
}

/** Main function */
int main(void)
{
    int n;
    n = 256;

    if (n <= 0 || (n & (n - 1)) != 0)
    {
        printf("Matrix size must be a positive power of 2!\n");
        return 1;
    }

    int **a = allocate_matrix(n);
    int **b = allocate_matrix(n);
    int **result = allocate_matrix(n);

    // Populate matrices with random values
    srand(42); // Seed for reproducibility
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            a[i][j] = rand() % 10; // Random numbers in range [0, 9]
            b[i][j] = rand() % 10;
        }

    // Perform Strassen's matrix multiplication
    double start_time = omp_get_wtime();
    strassen(a, b, result, n);
    double end_time = omp_get_wtime();

    printf("Time taken: %f seconds\n", end_time - start_time);

    // Free allocated memory
    free_matrix(a, n);
    free_matrix(b, n);
    free_matrix(result, n);

    return 0;
}
/** @} */
