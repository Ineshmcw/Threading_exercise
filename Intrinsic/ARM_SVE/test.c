#include <arm_sve.h>
#include <stdio.h>

void test_zip1(int32_t *input1, int32_t *input2, int32_t *output1, int32_t *output2) {
    svbool_t pg = svptrue_b32();
    svbool_t mask1 = svdupq_b32(1,1,0,0); // First two elements active
    // svbool_t mask2 = svdupq_b32(1,1,0,0);
    // svbool_t pred = svwhilelt_b32(2, 4); // Set first two false, last two true


    // Load vectors A and B
    svint32_t a = svld1(pg, input1);
    svint32_t b = svld1(pg, input2);

    // Interleave first two elements from each vector
    svint32_t result1 = svsplice(mask1, a, b);
    // svint32_t result2 = svsplice(mask2, a, b);

    // Store result
    svst1(pg, output1, result1);
    // svst1(pg, output2, result2);
}

int main() {
    int32_t input1[] = {1, 2, 3, 4}; // Vector A
    int32_t input2[] = {5, 6, 7, 8}; // Vector B
    int32_t output1[4] = {0}; // Output buffer
    int32_t output2[4] = {0}; // Output buffer

    // Call function
    test_zip1(input1, input2, output1, output2);

    // Print result
    printf("Result: [%d, %d, %d, %d]\n", output1[0], output1[1], output1[2], output1[3]);
    printf("Result: [%d, %d, %d, %d]\n", output2[0], output2[1], output2[2], output2[3]);


    return 0;
}




// #include <stdio.h>
// #include <arm_sve.h>

// void apply_sve_operation(float* matrix[4][4], float* result[4][4]) {
//     // Define masks
//     svbool_t pred = svptrue_b32();
//     svbool_t mask1 = svdupq_b32(1,1,0,0); // First two elements active
//     svbool_t mask2 = svdupq_b32(0,0,1,1); // First two elements active

//     // Load rows into vector registers
//     svfloat32_t v0 = svld1_f32(pred, matrix[0]);
//     svfloat32_t v1 = svld1_f32(pred, matrix[1]);
//     svfloat32_t v2 = svld1_f32(pred, matrix[2]);
//     svfloat32_t v3 = svld1_f32(pred, matrix[3]);

//     // Apply svsplice using mask1
//     svfloat32_t r0 = svsplice_f32(mask1, v0, v2);
//     svfloat32_t r1 = svsplice_f32(mask1, v1, v3);

//     // Apply svsplice using mask2
//     svfloat32_t r2 = svsplice_f32(mask2, v0, v2);
//     svfloat32_t r3 = svsplice_f32(mask2, v1, v3);

//     float res1[4] = {0};
//     float res2[4] = {0};
//     float res3[4] = {0};
//     float res4[4] = {0};
//     // Store the results back to the result matrix
//     svst1_f32(pred, res1, r0);
//     svst1_f32(pred, res2, r1);
//     svst1_f32(pred, res3, r2);
//     svst1_f32(pred, res4, r3);

//     // printf("Result matrix:\n");
//     // for (int i = 0; i < 4; i++) {
//     //     for (int j = 0; j < 4; j++) {
//     //         printf("%.1f ", result[i][j]);
//     //     }
//     //     printf("\n");
//     // }
//     printf("res1%.1f",res1);
//     printf("res1%.1f",res2);
//     printf("res1%.1f",res3);
//     printf("res1%.1f",res4);
// }

// int main() {
//     float matrix[4][4] = {
//         {1, 2, 3, 4},
//         {5, 6, 7, 8},
//         {9, 10, 11, 12},
//         {13, 14, 15, 16}
//     };

//     float result[4][4];

//     apply_sve_operation(matrix, result);

//     // printf("Result matrix:\n");
//     // for (int i = 0; i < 4; i++) {
//     //     for (int j = 0; j < 4; j++) {
//     //         printf("%.1f ", result[i][j]);
//     //     }
//     //     printf("\n");
//     // }

//     return 0;
// }
