/**
 * @file
 * @brief Implementation of [merge
 * sort](https://en.wikipedia.org/wiki/Merge_sort) algorithm
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

/**
 * @addtogroup sorting Sorting algorithms
 * @{
 */
/** Swap two integer variables
 * @param [in,out] a pointer to first variable
 * @param [in,out] b pointer to second variable
 */
void swap(int *a, int *b)
{
    int t;
    t = *a;
    *a = *b;
    *b = t;
}

/**
 * @brief Perform merge of segments.
 *
 * @param a array to sort
 * @param l left index for merge
 * @param r right index for merge
 * @param n total number of elements in the array
 */
void merge(int *a, int l, int r, int n)
{
    int *b = (int *)malloc(n * sizeof(int)); /* dynamic memory must be freed */
    if (b == NULL)
    {
        printf("Can't Malloc! Please try again.");
        exit(EXIT_FAILURE);
    }
    int c = l;
    int p1, p2;
    p1 = l;
    p2 = ((l + r) / 2) + 1;
    while ((p1 < ((l + r) / 2) + 1) && (p2 < r + 1))
    {
        if (a[p1] <= a[p2])
        {
            b[c++] = a[p1];
            p1++;
        }
        else
        {
            b[c++] = a[p2];
            p2++;
        }
    }

    if (p2 == r + 1)
    {
        while ((p1 < ((l + r) / 2) + 1))
        {
            b[c++] = a[p1];
            p1++;
        }
    }
    else
    {
        while ((p2 < r + 1))
        {
            b[c++] = a[p2];
            p2++;
        }
    }

    for (c = l; c < r + 1; c++) a[c] = b[c];

    free(b);
}

/** Merge sort algorithm implementation
 * @param a array to sort
 * @param n number of elements in the array
 * @param l index to sort from
 * @param r index to sort till
 */
void merge_sort(int *a, int n, int l, int r)
{
    if (r - l == 1)
    {
        if (a[l] > a[r])
            swap(&a[l], &a[r]);
    }
    else if (l != r)
    {
        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp task
                {
                    merge_sort(a, n, l, (l + r) / 2);
                }
                #pragma omp task
                {
                    merge_sort(a, n, ((l + r) / 2) + 1, r);
                }
            }
        }
        merge(a, l, r, n);
    }

    /* no change if l == r */
}
/** @} */

/** Main function */
int main(void)
{
    int *a, n, i;
    printf("Enter Array size: ");
    // scanf("%d", &n);
    n = 5000000;
    if (n <= 0) /* exit program if arraysize is not greater than 0 */
    {
        printf("Array size must be greater than 0!\n");
        return 1;
    }

    a = (int *)malloc(n * sizeof(int));
    if (a == NULL) /* exit program if can't malloc memory */
    {
        printf("Can't allocate memory! Please try again.\n");
        return 1;
    }

    /* Seed random number generator */
    srand(time(NULL));

    /* Populate the array with random integers */
    for (i = 0; i < n; i++)
    {
        a[i] = rand() % 1000; // Random numbers in range [0, 999]
    }

    // printf("Generated Array: ");
    // for (i = 0; i < n; i++)
    // {
    //     printf("%d ", a[i]);
    // }
    printf("\n");

    /* Call merge_sort */
    double start_time, end_time;
    start_time = omp_get_wtime();
    merge_sort(a, n, 0, n - 1);
    end_time = omp_get_wtime();

    printf("Time taken: %f\n", end_time - start_time);
    // printf("Sorted Array: ");
    // for (i = 0; i < n; i++)
    // {
    //     printf("%d ", a[i]);
    // }
    printf("\n");

    free(a);

    return 0;
}
