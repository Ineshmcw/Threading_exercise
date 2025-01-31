#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <math.h>
#include <random>

using namespace std;

const int N = 128;

// Standard matrix multiplication using 2D arrays
void mulMat(float m1[N][N], float m2[N][N], float res[N][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float sum = 0.0f;
            for (int k = 0; k < n; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            res[i][j] = sum;
        }
    }
}


// Function to check if two matrices are equal element-wise
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






bool areMatricesEqual(float m1[N][N], float m2[N][N], float tolerance = 1e-6) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Compare the absolute difference with tolerance
            float diff = fabs(m1[i][j] - m2[i][j]);
            if (diff > tolerance) {
                cout << "m1[" << i << "][" << j << "] = " << m1[i][j] << ", m2[" << i << "][" << j << "] = " << m2[i][j] << endl;
                cout << "Difference: " << diff << endl;
                return false; // Matrices are not equal
            }
        }
    }
    cout << "Matrices are equal within the tolerance." << endl;
    return true; // Matrices are equal within the tolerance
}

void printMatrices(float m1[N][N], float m2[N][N]) {
    cout << "Matrix 1:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << m1[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Matrix 2:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << m2[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    float m1[N][N], m2[N][N], res1[N][N] = {0}, res2[N][N] = {0};

    // Random number generator setup
    // random_device rd;
    // mt19937 gen(rd());
    // uniform_real_distribution<float> dis(0.0f, 10.0f); // Random values between 0 and 10

    // // Initialize matrices with random numbers
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         m1[i][j] = dis(gen); // Fill m1 with random numbers
    //         m2[i][j] = dis(gen); // Fill m2 with random numbers
    //     }
    // }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            m1[i][j] = i * N + j + 1;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            m2[i][j] = (N * N) - (i * N + j);
        }
    }

    // Measure time for non-AVX implementation using 2D arrays
    auto start = chrono::high_resolution_clock::now();
    mulMat(m1, m2, res1, N);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_non_avx = end - start;

    cout << "Execution time (Before AVX): " << elapsed_non_avx.count() << " seconds\n";

    // Measure time for AVX-based implementation
    start = chrono::high_resolution_clock::now();
    mulMatAVX(m1, m2, res2, N);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_avx = end - start;

    printMatrices(m1, m2);
    printMatrices(res1, res2);


    cout << "Execution time (With AVX): " << elapsed_avx.count() << " seconds\n";
    double speed_difference = ((elapsed_non_avx.count() - elapsed_avx.count()) / elapsed_non_avx.count()) * 100.0;
    cout << "Speed improvement with AVX: " << speed_difference << "%\n";

    // Now check if the results are equal
    if (areMatricesEqual(res1, res2)) {
        cout << "The results of the two matrix multiplication functions are equal!" << endl;
    } else {
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
 * This function performs matrix multiplication using AVX (Advanced Vector Extensions) intrinsics.
 * Unlike the FMA-based version, this implementation stores the intermediate sum results in a temporary array 
 * before accumulating them into a scalar. This approach results in a **77% speed improvement for N = 128**.
 * 
 * Optimizations used:
 * 1. **Vectorized Loading:** Loads 8 floating-point values at a time using `_mm256_loadu_ps`.
 * 2. **Vectorized Multiplication and Summation:** Uses `_mm256_mul_ps` for element-wise multiplication and `_mm256_add_ps` to accumulate results.
 * 3. **Temporary Storage:** Stores the AVX register result into an array using `_mm256_storeu_ps`, then sums up the 8 elements sequentially.
 * 4. **Matrix Transposition:** The second matrix is transposed using `unpacklo`, `unpackhi`, `shuffle`, and `permute` instructions to improve memory access patterns.
 * 
 * Performance:
 * - Transposing the matrix with '__m256_loadu_ps' and storing column wise manually, results in **69% improvement for N = 128**.
 * - Achieves approximately **77% improvement for N = 128** ( with shuffle and permute transpose).
 * - While slightly slower than the `_mm256_fmadd_ps` version, this method remains efficient.
 * 
 * Parameters:
 * - `m1`: Input matrix of size N x N.
 * - `trans`: Transposed second matrix of size N x N (precomputed for better cache efficiency).
 * - `res`: Output matrix of size N x N to store the result.
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
 * This function performs matrix multiplication using AVX (Advanced Vector Extensions) intrinsics for optimized performance. 
 * The input matrices are assumed to be of size N x N, where N is a globally defined constant and a multiple of 8.
 * The function utilizes AVX2 instructions for efficient vectorized computation.
 * 
 * Optimizations used:
 * 1. **Vectorized Loading:** Loads 8 floating-point values at a time using `_mm256_loadu_ps`.
 * 2. **Fused Multiply-Add (FMA):** Uses `_mm256_fmadd_ps` to perform `sum += vecA * vecB` in a single instruction, reducing instruction count and improving efficiency.
 * 3. **Matrix Transposition:** The second matrix is transposed before multiplication to improve cache locality and avoid non-contiguous memory access.
 * 4. **Horizontal Summation:** Uses `_mm256_hadd_ps` to efficiently sum partial results across vector lanes.
 * 5. **Final Summation and Extraction:** Converts the 256-bit AVX register into two 128-bit registers, sums them up, and extracts the final scalar value.
 * 
 * Performance Improvement:
 * - Transposing the matrix with '__m256_loadu_ps' and storing column wise manually, the performance improvement is around **74% for N = 128**.
 * - Using `_mm256_loadu_ps`, `_mm256_fmadd_ps`, `_mm256_hadd_ps`, and transposition, the performance improvement is approximately **81% for N = 128**.
 * 
 * Parameters:
 * - `m1`: Input matrix of size N x N.
 * - `trans`: Transposed second matrix of size N x N (precomputed for better cache efficiency).
 * - `res`: Output matrix of size N x N to store the result.
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
 * - **74% speedup for N = 128** using transposition with manual column-wise storage.
 * - **81% speedup for N = 128** using transposition with unpacklo and unpackhi.
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
