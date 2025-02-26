#include <stdio.h>
#include <arm_neon.h>


// To run
// arm-linux-gnueabihf-gcc -march=armv7-a -mfpu=neon -mfloat-abi=hard -o exercise exercise.c
// qemu-arm -L /usr/arm-linux-gnueabihf/ ./exercise

float dot_product(float *A, float *B, int N) {
    float32x4_t sum = vmovq_n_f32(0.0);
    for (int i = 0; i < N; i+=4) {
        float32x4_t a = vld1q_f32(A+i);
        float32x4_t b = vld1q_f32(B+i);
        float32x4_t c_i = vmulq_f32(a,b);
        sum = vaddq_f32(sum, c_i);
    }
    float32x2_t sum_low, sum_high, final_sum;

    sum_low = vget_low_f32(sum);
    sum_high = vget_high_f32(sum);
    final_sum = vadd_f32(sum_low, sum_high);

    return vget_lane_f32(vpadd_f32(final_sum, final_sum), 0); 
}

int main() {
    float A[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    float B[8] = {8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

    float result = dot_product(A, B, 8);

    printf("Dot Product Result: %f\n", result);

    return 0;
}


/* Methord 1 */
// float dot_product(float *A, float *B, int N) {
//     float32x4_t sum = vmovq_n_f32(0.0);
//     for (int i = 0; i < N; i+=4) {
//         float32x4_t a = vld1q_f32(A+i);
//         float32x4_t b = vld1q_f32(B+i);
//         float32x4_t c_i = vmulq_f32(a,b);
//         sum = vaddq_f32(sum, c_i);
//     }
//     float final_sum[4];
//     vst1q_f32(final_sum, sum);
//     return final_sum[0] + final_sum[1] + final_sum[2] + final_sum[3];
// }


/* Methord 2 */
// float dot_product(float *A, float *B, int N) {
//     float32x4_t sum = vmovq_n_f32(0.0);
//     for (int i = 0; i < N; i+=4) {
//         float32x4_t a = vld1q_f32(A+i);
//         float32x4_t b = vld1q_f32(B+i);
//         float32x4_t c_i = vmulq_f32(a,b);
//         sum = vaddq_f32(sum, c_i);
//     }
//     float32x2_t sum_low, sum_high, final_sum;
//     sum_low = vget_low_f32(sum);
//     sum_high = vget_high_f32(sum);
//     final_sum = vadd_f32(sum_low, sum_high);
//     return vget_lane_f32(vpadd_f32(final_sum, final_sum), 0); 
// }


/* Multiplication */
// void vector_multiplication(float *A, float *B, float *C, int N) {
//     for (int i = 0; i < N; i+=4) {
//         float32x4_t a = vld1q_f32(A+i);
//         float32x4_t b = vld1q_f32(B+i);
//         float32x4_t c_i = vmulq_f32(a,b);
//         vst1q_f32(C + i,c_i);
//     }
// }


/* Addition */
// void vector_addition(float *A, float *B, float *C, int N) {


//     for (int i = 0; i < N; i += 4) {
//         float32x4_t a = vld1q_f32(A + i);
//         float32x4_t b = vld1q_f32(B + i);
//         float32x4_t c_i = vaddq_f32(a, b);
        
//         vst1q_f32(C + i,c_i);
//     }
// }