#include <stdio.h>
#include <arm_neon.h>

void print_f32x4(float32x4_t vec, const char* label) {
    float vals[4];
    vst1q_f32(vals, vec);  // Store the vector into an array
    printf("%s: { %.1f, %.1f, %.1f, %.1f }\n", label, vals[0], vals[1], vals[2], vals[3]);
}

int main() {
    // Define two input vectors
    float32x4_t a = {1.0, 2.0, 3.0, 4.0};
    float32x4_t b = {5.0, 6.0, 7.0, 8.0};

    // Perform interleaving: extract the second half
    float32x4_t result = vget_low_f32(a, b);

    // Print results
    print_f32x4(a, "a");
    print_f32x4(b, "b");
    print_f32x4(result, "vzip2q_f32(a, b)");

    return 0;
}
