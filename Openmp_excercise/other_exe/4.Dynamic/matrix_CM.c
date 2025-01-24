/**
 * @file
 * @brief [Matrix Chain
 * Order](https://en.wikipedia.org/wiki/Matrix_chain_multiplication)
 * @details
 * From Wikipedia: Matrix chain multiplication (or the matrix chain ordering
 * problem) is an optimization problem concerning the most efficient way to
 * multiply a given sequence of matrices. The problem is not actually to perform
 * the multiplications, but merely to decide the sequence of the matrix
 * multiplications involved.
 * @author [CascadingCascade](https://github.com/CascadingCascade)
 */

#include <assert.h>  /// for assert
#include <limits.h>  /// for INT_MAX macro
#include <stdio.h>   /// for IO operations
#include <stdlib.h>  /// for malloc() and free()
#include <omp.h>

/**
 * @brief Finds the optimal sequence using the classic O(n^3) algorithm.
 * @param l length of cost array
 * @param p costs of each matrix
 * @param s location to store results
 * @returns number of operations
 */
int matrixChainOrder(int l, const int *p, int *s)
{
    // mat stores the cost for a chain that starts at i and ends on j (inclusive
    // on both ends)
    int **mat = malloc(l * sizeof(int *));
    for (int i = 0; i < l; ++i)
    {
        mat[i] = malloc(l * sizeof(int));
    }

    for (int i = 0; i < l; ++i)
    {
        mat[i][i] = 0;
    }
    // cl denotes the difference between start / end indices, cl + 1 would be
    // chain length.
    for (int cl = 1; cl < l; ++cl)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < l - cl; ++i)
        {

            int j = i + cl;
            mat[i][j] = INT_MAX;
            
            for (int div = i; div < j; ++div)
            {
                int q = mat[i][div] + mat[div + 1][j] + p[i] * p[div] * p[j];
                if (q < mat[i][j])
                {
                    mat[i][j] = q;
                    s[i * l + j] = div;
                }
            
            }
        }
        
        
    }
    
    
    
    int result = mat[0][l - 1];

    // Free dynamically allocated memory
    for (int i = 0; i < l; ++i)
    {
        free(mat[i]);
    }
    free(mat);

    return result;
}

/**
 * @brief Recursively prints the solution
 * @param l dimension of the solutions array
 * @param s solutions
 * @param i starting index
 * @param j ending index
 * @returns void
 */
void printSolution(int l, int *s, int i, int j)
{
    if (i == j)
    {
        printf("A%d", i);
        return;
    }
    putchar('(');
    printSolution(l, s, i, s[i * l + j]);
    printSolution(l, s, s[i * l + j] + 1, j);
    putchar(')');
}

int main() {
    int len = 1000;
    int *sizes = malloc(len * sizeof(int));
    int *sol = malloc(len * len * sizeof(int));

    for (int i = 0; i < len; ++i) {
        sizes[i] = rand() % 100 + 1;
    }

    double start = omp_get_wtime();
    int result = matrixChainOrder(len, sizes, sol);
    double end = omp_get_wtime();

    printf("Result: %d", result);
    // printf("Optimal ordering: ");
    // printSolution(len, sol, 0, len - 1);
    printf("\nTime taken: %f seconds\n", end - start);

    free(sizes);
    free(sol);

    return 0;
}