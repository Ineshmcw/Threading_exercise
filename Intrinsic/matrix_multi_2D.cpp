#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <math.h>
#include <random>
#include <stdio.h>

using namespace std;

const int N = 123;
const int M = 124;

// Standard matrix multiplication using 2D arrays
// void mulMat(float m1[N][N], float m2[N][N], float res[N][N], int n) {
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

void _print_m256(__m256 val)
{
    float *f = (float *)&val;
    printf("{%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f}\n",
           f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
}

void printMatrices(float res[M][M], bool isAVX = false)
{
    if (isAVX)
        cout << "Result Matrix(after AVX):" << endl;
    else
        cout << "Result Matrix(before AVX):" << endl;
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            cout << res[i][j] << " ";
        }
        cout << endl;
    }
}

void mulMatAVX(float m1[M][N], float m2[N][M], float res[M][M], int m, int n)
{

    int chunk_size = 8;
    int aligned_n = (n / chunk_size) * chunk_size; // Largest multiple of 8 ≤ N
    int aligned_m = (m / chunk_size) * chunk_size; // Largest multiple of 8 ≤ N

    // Process full 8×8 blocks
    for (int bi = 0; bi < aligned_m; bi += chunk_size)
    {
        for (int bj = 0; bj < aligned_m; bj += chunk_size)
        {
            for (int bk = 0; bk < aligned_n; bk += chunk_size)
            {

                // Load 8 rows of m2
                for (int rowm2 = bk; rowm2 < bk + chunk_size; rowm2 += 8)
                {
                    __m256 row0 = _mm256_loadu_ps(&m2[rowm2 + 0][bj]);
                    __m256 row1 = _mm256_loadu_ps(&m2[rowm2 + 1][bj]);
                    __m256 row2 = _mm256_loadu_ps(&m2[rowm2 + 2][bj]);
                    __m256 row3 = _mm256_loadu_ps(&m2[rowm2 + 3][bj]);
                    __m256 row4 = _mm256_loadu_ps(&m2[rowm2 + 4][bj]);
                    __m256 row5 = _mm256_loadu_ps(&m2[rowm2 + 5][bj]);
                    __m256 row6 = _mm256_loadu_ps(&m2[rowm2 + 6][bj]);
                    __m256 row7 = _mm256_loadu_ps(&m2[rowm2 + 7][bj]);

                    for (int rowm1 = bi; rowm1 < bi + chunk_size; rowm1++)
                    {
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
            // load the remaning rows of m2
            for (int rowm2 = aligned_n; rowm2 < n; rowm2++)
            {
                __m256 row = _mm256_loadu_ps(&m2[rowm2][bj]);

                for (int rowm1 = bi; rowm1 < bi + chunk_size; rowm1++)
                {
                    __m256 res_store = _mm256_loadu_ps(&res[rowm1][bj]);
                    res_store = _mm256_add_ps(res_store, _mm256_mul_ps(_mm256_set1_ps(m1[rowm1][rowm2]), row));
                    _mm256_storeu_ps(&res[rowm1][bj], res_store);
                }
            }
        }
    }

    // Handling remaining rows and columns
    for (int i = 0; i < m; i++)
    {
        for (int j = aligned_n; j < m; j++)
        { // Remaining columns
            float sum = 0;
            for (int k = 0; k < n; k++)
            {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }

    for (int i = aligned_m; i < m; i++)
    { // Remaining rows
        for (int j = 0; j < aligned_n; j++)
        {
            float sum = 0;
            for (int k = 0; k < n; k++)
            {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
}

bool areMatricesEqual(float m1[M][M], float m2[M][M], float tolerance = 1e-6)
{
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; ++j)
        { // Use M instead of N here
            float diff = fabs(m1[i][j] - m2[i][j]);
            if (diff > tolerance)
            {
                cout << "m1[" << i << "][" << j << "] = " << m1[i][j] << ", m2[" << i << "][" << j << "] = " << m2[i][j] << endl;
                cout << "Difference: " << diff << endl;
                return false; // Matrices are not equal
            }
        }
    }
    cout << "Matrices are equal within the tolerance." << endl;
    return true; // Matrices are equal within the tolerance
}

int main()
{
    // float m1[N][N], m2[N][N], res1[N][N] = {0}, res2[N][N] = {0};
    float m1[M][N] __attribute__((aligned(32)));
    float m2[N][M] __attribute__((aligned(32)));

    float res1[M][M] = {0}, res2[M][M] = {0};

    cout << "m1 matrix:" << endl;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            m1[i][j] = 1;
            cout << m1[i][j] << " ";
        }
        cout << endl;
    }

    cout << "m2 matrix:" << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            m2[i][j] = 2;
            cout << m2[i][j] << " ";
        }
        cout << endl;
    }

    // Measure time for non-AVX implementation using 2D arrays
    auto start = chrono::high_resolution_clock::now();
    mulMat(m1, m2, res1, M, N);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_non_avx = end - start;

    // Measure time for AVX-based implementation
    start = chrono::high_resolution_clock::now();
    mulMatAVX(m1, m2, res2, M, N);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_avx = end - start;

    printMatrices(res1);
    printMatrices(res2, true);

    cout << "Execution time (Before AVX): " << elapsed_non_avx.count() << " seconds\n";
    cout << "Execution time (With AVX): " << elapsed_avx.count() << " seconds\n";
    double speed_difference = ((elapsed_non_avx.count() - elapsed_avx.count()) / elapsed_non_avx.count()) * 100.0;
    cout << "Speed improvement with AVX: " << speed_difference << "%\n";

    // Now check if the results are equal
    if (areMatricesEqual(res1, res2))
    {
        cout << "The results of the two matrix multiplication functions are equal!" << endl;
    }
    else
    {
        cout << "The results of the two matrix multiplication functions are NOT equal!" << endl;
    }

    return 0;
}

/*

Transposing matrix with loadu and manually storing column wise.

*/

/*
    for (int rows = 0; rows < N; rows++) {
        for (int columns = 0; columns < N; columns += 8) {
            __m256 row = _mm256_loadu_ps(&m2[rows][columns]);
            for (int k = 0; k < 8; k++) {
                trans[columns + k][rows] = ((float*)&row)[k];
            }
        }
    }
*/

/*
Transposing with unpacklo, unpackhi, shuffle and permute
*/

/*

    float trans[N][N] = {0};
    for (int i = 0; i < N; i += 8) {
        for (int j = 0; j < N; j += 8) {
                __m256 row0 = _mm256_loadu_ps(&m2[i + 0][j]);
                __m256 row1 = _mm256_loadu_ps(&m2[i + 1][j]);
                __m256 row2 = _mm256_loadu_ps(&m2[i + 2][j]);
                __m256 row3 = _mm256_loadu_ps(&m2[i + 3][j]);
                __m256 row4 = _mm256_loadu_ps(&m2[i + 4][j]);
                __m256 row5 = _mm256_loadu_ps(&m2[i + 5][j]);
                __m256 row6 = _mm256_loadu_ps(&m2[i + 6][j]);
                __m256 row7 = _mm256_loadu_ps(&m2[i + 7][j]);

                // Transpose using AVX intrinsics
                __m256 t0 = _mm256_unpacklo_ps(row0, row1);
                __m256 t1 = _mm256_unpackhi_ps(row0, row1);
                __m256 t2 = _mm256_unpacklo_ps(row2, row3);
                __m256 t3 = _mm256_unpackhi_ps(row2, row3);
                __m256 t4 = _mm256_unpacklo_ps(row4, row5);
                __m256 t5 = _mm256_unpackhi_ps(row4, row5);
                __m256 t6 = _mm256_unpacklo_ps(row6, row7);
                __m256 t7 = _mm256_unpackhi_ps(row6, row7);

                __m256 tt0 = _mm256_shuffle_ps(t0, t2, 0x44);
                __m256 tt1 = _mm256_shuffle_ps(t0, t2, 0xEE);
                __m256 tt2 = _mm256_shuffle_ps(t1, t3, 0x44);
                __m256 tt3 = _mm256_shuffle_ps(t1, t3, 0xEE);
                __m256 tt4 = _mm256_shuffle_ps(t4, t6, 0x44);
                __m256 tt5 = _mm256_shuffle_ps(t4, t6, 0xEE);
                __m256 tt6 = _mm256_shuffle_ps(t5, t7, 0x44);
                __m256 tt7 = _mm256_shuffle_ps(t5, t7, 0xEE);

                __m256 r0 = _mm256_permute2f128_ps(tt0, tt4, 0x20);
                __m256 r1 = _mm256_permute2f128_ps(tt1, tt5, 0x20);
                __m256 r2 = _mm256_permute2f128_ps(tt2, tt6, 0x20);
                __m256 r3 = _mm256_permute2f128_ps(tt3, tt7, 0x20);
                __m256 r4 = _mm256_permute2f128_ps(tt0, tt4, 0x31);
                __m256 r5 = _mm256_permute2f128_ps(tt1, tt5, 0x31);
                __m256 r6 = _mm256_permute2f128_ps(tt2, tt6, 0x31);
                __m256 r7 = _mm256_permute2f128_ps(tt3, tt7, 0x31);

                // Store transposed data
                _mm256_storeu_ps(&trans[j + 0][i], r0);
                _mm256_storeu_ps(&trans[j + 1][i], r1);
                _mm256_storeu_ps(&trans[j + 2][i], r2);
                _mm256_storeu_ps(&trans[j + 3][i], r3);
                _mm256_storeu_ps(&trans[j + 4][i], r4);
                _mm256_storeu_ps(&trans[j + 5][i], r5);
                _mm256_storeu_ps(&trans[j + 6][i], r6);
                _mm256_storeu_ps(&trans[j + 7][i], r7);
            }
        }

*/

/*

Used loadu and storeu and used temp variable to store the results of the multiplication
Along with transposing with unpacklo, unpackhi, shuffle and permute speed imporovement is 77% for N=128

*/

/*
 * Function: matrix_multiply_avx_temp
 *
 * Performs matrix multiplication using AVX intrinsics, storing intermediate sums in a temporary array.
 *
 * Optimizations:
 * 1. **Vectorized Loading & Computation:** Uses `_mm256_loadu_ps`, `_mm256_mul_ps`, and `_mm256_add_ps`.
 * 2. **Temporary Storage:** Stores results with `_mm256_storeu_ps` before summing.
 * 3. **Matrix Transposition:** Improves cache efficiency using shuffle and permute operations.
 *
 * Performance:
 * - **69% improvement for N = 128** with manual transposition.
 * - **77% improvement for N = 128** using shuffle and permute.
 *
 * Parameters:
 * - `m1`: Input matrix (N x N).
 * - `trans`: Transposed second matrix (N x N).
 * - `res`: Output matrix (N x N).
 */

/*
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            __m256 sum = _mm256_setzero_ps();

            for (int k = 0; k < N; k += 8) {
                __m256 vecA = _mm256_loadu_ps(&m1[i][k]);
                __m256 vecB = _mm256_loadu_ps(&trans[j][k]);
                __m256 prod = _mm256_mul_ps(vecA, vecB);
                sum = _mm256_add_ps(sum, prod);
            }

            float temp[8];
            _mm256_storeu_ps(temp, sum);
            res[i][j] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
        }
    }

*/

/*

Used loadu, fmadd, hadd(horizontal add) and used castps, extract(sum the low and high)
Along with transposing with unpacklo, unpackhi, shuffle and permute speed imporovement is 81% for N=128

*/

/*
 * Function: matrix_multiply_avx
 *
 * Performs N x N matrix multiplication using AVX2 intrinsics for optimized performance.
 *
 * Optimizations:
 * 1. **Vectorized Loading:** Uses `_mm256_loadu_ps` to load 8 values at a time.
 * 2. **FMA Operations:** Uses `_mm256_fmadd_ps` for efficient fused multiply-add computations.
 * 3. **Matrix Transposition:** Pre-transposes `m2` to improve cache locality.
 * 4. **Efficient Summation:** Utilizes `_mm256_hadd_ps` for fast horizontal summation.
 *
 * Performance:
 * - **74% speedup for N = 128** with manual transposition.
 * - **81% speedup for N = 128** using using shuffle and permute.
 *
 * Parameters:
 * - `m1`: Input matrix (N x N).
 * - `trans`: Pre-transposed second matrix (N x N).
 * - `res`: Output matrix (N x N).
 */

/*
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            __m256 sum = _mm256_setzero_ps();

            for (int k = 0; k < N; k += 8) {
                __m256 vecA = _mm256_loadu_ps(&m1[i][k]);
                __m256 vecB = _mm256_loadu_ps(&trans[j][k]); // transposed matrix ( m2 )

                sum = _mm256_fmadd_ps(vecA, vecB, sum);
            }

            sum = _mm256_hadd_ps(sum, sum);
            sum = _mm256_hadd_ps(sum, sum);
            __m128 sum_low = _mm256_castps256_ps128(sum);
            __m128 sum_high = _mm256_extractf128_ps(sum, 1);
            __m128 final_sum = _mm_add_ps(sum_low, sum_high);
            res[i][j] = _mm_cvtss_f32(final_sum);
        }
    }
*/

/*
 * Function: mulMatAVX
 *
 * Performs matrix multiplication using AVX intrinsics for efficient computation.
 *
 * Optimizations:
 * 1. **Vectorized Column Extraction:** Stores `m2` columns in a temporary array to improve cache locality.
 * 2. **AVX Operations:** Uses `_mm256_loadu_ps`, `_mm256_mul_ps`, and `_mm256_add_ps` for efficient SIMD computation.
 * 3. **Efficient Summation:** Stores intermediate sums in a temporary array before accumulating results.
 *
 * Performance:
 * - **65% speedup for N = 128** using AVX.
 *
 * Parameters:
 * - `m1`, `m2`: Input matrices of size N x N.
 * - `res`: Output matrix of size N x N.
 * - `n`: Matrix size (assumed to be a multiple of 8 for AVX).
 */

/*

void mulMatAVX(float m1[N][N], float m2[N][N], float res[N][N], int n) {
    for (int columns = 0; columns < n; columns++) {

        float column[N];
        for (size_t row = 0; row < N; row ++ ) {
            column[row] = m2[row][columns];
        }

        for (int row = 0; row < N; row++) {
            __m256 sum = _mm256_setzero_ps();

            for (int k = 0; k < N; k += 8) {
                __m256 vecA = _mm256_loadu_ps(&m1[row][k]);
                __m256 vecB = _mm256_loadu_ps(&column[k]);
                __m256 prod = _mm256_mul_ps(vecA, vecB);
                sum = _mm256_add_ps(sum, prod);
            }

            float temp[8];
            _mm256_storeu_ps(temp, sum);
            res[row][columns] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
        }
    }
}

*/

/*
New methord logic
*/

/*
    //1. For nxn for 256 intrinsic, first we load the 8 rows of matrix 2.
    //2. we take in the first element of matrix 1 and muliply with first row of matrix 2
    //3. we store the each multiplied values in the res matrix first row respectively.
    //4. we then modifiy step 2 by taking the 2nd element of 1st matrix 1st row and multiply them with the second row of matrix 2.
    //5. store them in by adding the res matrix first row previously saved values.
    //6. we continue this step until all the elements of the first row is properly identified.
    //7. now we move to the second row of first matrix and proceed the same.
*/

/*
Starting idea for block wise matrix multiplication.
*/

/*

    for(int i = 0; i < n; i += 8) {
        for(int j = 0; j < n; j += 8) {
            for(int k = 0; k < n; k += 8) {
                for(int x = i; x < i + 8; x++) {
                    for(int y = j; y < j + 8; y++) {
                        __m256 c_vec = _mm256_setzero_ps();
                        for(int z = k; z < k + 8; z += 8) {
                            __m256 a_vec = _mm256_loadu_ps(&m1[i][k]);  // Load 8 elements of A[i][k]
                            __m256 b_vec = _mm256_loadu_ps(&m2[k][j]);  // Load 8 elements of B[k][j]
                            c_vec = _mm256_fmadd_ps(a_vec, b_vec, c_vec); // FMA: c_vec += a_vec * b_vec
                        }
                        res[i + x][j + y] = sum;
                    }
                }
            }
        }
    }
*/

// for(int rows = 0; rows < n; rows += 50) {
//     for(int cols = 0; cols < n; cols += 8) {
//         for(int col8 = cols; col8 < cols + 8; col8++) {
//             _m256 row1 = _mm256_loadu_ps(&m2[rows][col8]);
//             for(int row8 = rows; row8 < rows + 8; row8++) {

//         __m256

// // Function to check if two matrices are equal element-wise
// void mulMatAVX(float m1[N][N], float m2[N][N], float res[N][N], int n) {
//     for (int columns = 0; columns < n; columns++) { // traverse through each column of m2

//         float column[N] = {0};
//         for (size_t i = 0; i < N; i ++ ) {
//             _mm256_store_ps(&column[i], _mm256_loadu_ps(&m2[i][columns]));
//         }

//         for (int rows = 0; rows < n; rows++) { // traverse through each row of m1

//             __m256 sum = _mm256_setzero_ps();

//             for (int k = 0; k < n; k += 8) { //in 256 traverse through 8 elements at the same time and multiply them
//                 if (k + 7 < N)
//                 {
//                     __m256 data = _mm256_loadu_ps(&m1[rows][k]);
//                     __m256 data2 = _mm256_loadu_ps(&column[k]);
//                     sum = _mm256_add_ps(sum, _mm256_mul_ps(data, data2));
//                 }
//                 else
//                     for (int j = k; j < N; j++) {
//                         sum = _mm256_add_ps(sum, _mm256_set1_ps(m1[rows][j] * column[j]));
//                     }
//             }
//             __m256 total_sum = _mm256_setzero_ps();
//             __m256 temp;
//             _mm256_storeu_ps((float*)&temp, sum);
//             total_sum = _mm256_add_ps(total_sum, temp);

//             _mm256_storeu_ps((float*)&res[rows][columns], total_sum);
//         }
//     }
// }

// void mulMatAVX(float m1[N][N], float m2[N][N], float res[N][N]) {
//     for (int i = 0; i < N; i++) {   // Iterate over rows of A
//         for (int j = 0; j < N; j += 8) {  // Iterate over 8-element blocks
//             __m256 rowA1 = _mm256_loadu_ps(&m1[i][j]);
//             __m256 rowA2 = _mm256_loadu_ps(&m1[i][j+8]);

//             __m256 accumulator = _mm256_setzero_ps();  // Initialize the accumulator

//             for (int k = 0; k < N; k++) {  // Iterate over rows of B
//                 __m256 rowB1 = _mm256_loadu_ps(&m2[k][j]);
//                 __m256 rowB2 = _mm256_loadu_ps(&m2[k][j+8]);

//                 // Step 4: Blend A with zeros
//                 __m256 blendA1 = _mm256_blend_ps(rowA1, _mm256_setzero_ps(), 0b11110000);
//                 __m256 blendA2 = _mm256_blend_ps(rowA1, _mm256_setzero_ps(), 0b01000001);

//                 // Step 9: Blend next block of A
//                 __m256 blendA3 = _mm256_blend_ps(rowA2, _mm256_setzero_ps(), 0b11110000);
//                 __m256 blendA4 = _mm256_blend_ps(rowA2, _mm256_setzero_ps(), 0b01000001);

//                 // Step 5 & 10: Multiply blended A with B
//                 __m256 result1 = _mm256_mul_ps(blendA1, rowB1);
//                 __m256 result2 = _mm256_mul_ps(blendA3, rowB2);

//                 // Step 6 & 11: Accumulate results
//                 accumulator = _mm256_add_ps(accumulator, result1);
//                 accumulator = _mm256_add_ps(accumulator, result2);
//             }

//             // Step 16: Store the first 8 results of C
//             _mm256_storeu_ps(&res[i][j], accumulator);
//             // Step 17: Store the next 8 results of C
//             _mm256_storeu_ps(&res[i][j+8], accumulator);
//         }

//         // Step 18: Load the next row of A for the next iteration
//         __m256 rowA1_next = _mm256_loadu_ps(&m1[i+1][0]);
//         __m256 rowA2_next = _mm256_loadu_ps(&m1[i+1][8]);

//         // Step 19: Repeat steps 3-13 for the next row of B
//         for (int j = 0; j < N; j += 8) {
//             __m256 rowB1_next = _mm256_loadu_ps(&m2[i+1][j]);
//             __m256 rowB2_next = _mm256_loadu_ps(&m2[i+1][j+8]);

//             __m256 accumulator_next = _mm256_setzero_ps();  // Reset accumulator for next row
//             __m256 blendA1_next = _mm256_blend_ps(rowA1_next, _mm256_setzero_ps(), 0b11110000);
//             __m256 blendA3_next = _mm256_blend_ps(rowA2_next, _mm256_setzero_ps(), 0b11110000);

//             __m256 result1_next = _mm256_mul_ps(blendA1_next, rowB1_next);
//             __m256 result2_next = _mm256_mul_ps(blendA3_next, rowB2_next);

//             accumulator_next = _mm256_add_ps(accumulator_next, result1_next);
//             accumulator_next = _mm256_add_ps(accumulator_next, result2_next);

//             // Step 16: Store the next set of results in C
//             _mm256_storeu_ps(&res[i+1][j], accumulator_next);
//             _mm256_storeu_ps(&res[i+1][j+8], accumulator_next);
//         }
//     }
// }
