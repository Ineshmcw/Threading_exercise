#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

/*Swap function to swap two values*/
void swap(int *first, int *second)
{
    int temp = *first;
    *first = *second;
    *second = temp;
}

/*Partition method which selects a pivot
  and places each element which is less than the pivot value to its left
  and the elements greater than the pivot value to its right
  arr[] --- array to be partitioned
  lower --- lower index
  upper --- upper index
*/
int partition(int arr[], int lower, int upper)
{
    int i = (lower - 1);

    int pivot = arr[upper]; // Selects last element as the pivot value

    int j;
    for (j = lower; j < upper; j++)
    {
        if (arr[j] <= pivot)
        { // if current element is smaller than the pivot

            i++; // increment the index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }


swap(&arr[i + 1], &arr[upper]); // places the last element i.e, the pivot
                                // to its correct position

return (i + 1);
}

/*This is where the sorting of the array takes place
    arr[] --- Array to be sorted
    lower --- Starting index
    upper --- Ending index
*/
void quickSort(int arr[], int lower, int upper)
{
    if (upper > lower)
    {
        // partitioning index is returned by the partition method , partition
        // element is at its correct poition

        int partitionIndex = partition(arr, lower, upper);

        // Sorting elements before and after the partition index
        #pragma omp task shared(arr)
            quickSort(arr, lower, partitionIndex - 1);
        #pragma omp task shared(arr)
            quickSort(arr, partitionIndex + 1, upper);
    }
}

int main()
{
    int n,*arr,i;
    n = 1000000;

    arr = (int *)malloc(n * sizeof(int));
    if (arr == NULL) /* exit program if can't malloc memory */
    {
        printf("Can't allocate memory! Please try again.\n");
        return 1;
    }

    /* Seed random number generator */
    srand(time(NULL));

    /* Populate the array with random integers */
    for (i = 0; i < n; i++)
    {
        arr[i] = rand() % 1000; // Random numbers in range [0, 999]
    }
    double start_time, end_time;
    start_time = omp_get_wtime();
    #pragma omp parallel 
    {
        #pragma omp single 
            quickSort(arr, 0, n - 1);
    }
    end_time = omp_get_wtime();
    printf("Time taken: %f\n", end_time - start_time);
    free(arr);
    return 0;
}
