/*
 * To run:
 * aarch64-linux-gnu-gcc -g -O2 -march=armv8-a+sve -o exe exe.c
 * qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./exe
 */

#include <stdio.h>
#include <arm_sve.h>


const char* SkipWhiteSpace(const char* p, const char* end) {
    while (p != end && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')) {
    p++;
    }
    return p;
}

void exe(int *out, int *a, int *b, int *c, int N) {
    svbool_t predicate;
    for (int i = 0; i < N; i+= svcntw()) {
        predicate = svwhilelt_b32(i, N);

        svint32_t vec_a = svld1_s32(predicate, &a[i]);
        svint32_t vec_b = svld1_s32(predicate, &b[i]);
        svint32_t vec_c = svld1_s32(predicate, &c[i]);
        
        svst1_s32(predicate, &out[i], svmul_s32_m(predicate, svadd_s32_m( predicate, vec_a, vec_b), vec_c));
    }
}

int main() {
    int a[] = {1, 2, 3, 4};
    int b[] = {5, 6, 7, 8};
    int c[] = {9, 10, 11, 12};
    int N = sizeof(a) / sizeof(a[0]);
    int out[N];

    exe(out, a, b, c, N);

    for (int i = 0; i < N; i++) {
        printf("%d ", out[i]);
    }
    printf("\n");

    return 0;
}


/*

const char* SkipWhiteSpace(const char* p, const char* end) {
    while (p < end) {
        svbool_t pg = svwhilelt_b8((uint64_t)(p - end), (uint64_t)0);
        svuint8_t vec = svld1(pg, (const uint8_t*)p);

        // Compare against whitespace characters (space, '\n', '\r', '\t')
        svbool_t mask_space  = svcmpeq_n_u8(pg, vec, ' ');
        svbool_t mask_newline = svcmpeq_n_u8(pg, vec, '\n');
        svbool_t mask_carriage = svcmpeq_n_u8(pg, vec, '\r');
        svbool_t mask_tab = svcmpeq_n_u8(pg, vec, '\t');

        svbool_t mask_whitespace = svorr_b_z(pg, mask_space, mask_newline);
        mask_whitespace = svorr_b_z(pg, mask_whitespace, mask_carriage);
        mask_whitespace = svorr_b_z(pg, mask_whitespace, mask_tab);

        // Find first non-whitespace character
        int first_non_whitespace = svcntp_b8(svptrue_b8(), mask_whitespace);
        if (first_non_whitespace < svcntb()) {
            return p + first_non_whitespace;
        }

        // Move to the next chunk
        p += svcntb();
    }
    return p;
}


int main() {
    const char* text = "   \n\t  Hello, World!";
    const char* end = text + strlen(text); // Get the end of the string
    printf("character: %s\n",text);
    const char* result = SkipWhiteSpace(text, end);
    
    printf("First non-whitespace character: %s\n",result);
    printf("Remaining string: %s\n",result);

    return 0;
}

*/
// void exe(int *out, int *a, int *b, int *c, int N) {
//     svbool_t predicate;
//     for (int i = 0; i < N; i+= svcntw()) {
//         predicate = svwhilelt_b32(i, N);

//         svint32_t vec_a = svld1_s32(predicate, &a[i]);
//         svint32_t vec_b = svld1_s32(predicate, &b[i]);
//         svint32_t vec_c = svld1_s32(predicate, &c[i]);
        
//         svst1_s32(predicate, &out[i], svmul_s32_m(predicate, svadd_s32_m( predicate, vec_a, vec_b), vec_c));
//     }
// }

// int main() {
//     int a[] = {1, 2, 3, 4};
//     int b[] = {5, 6, 7, 8};
//     int c[] = {9, 10, 11, 12};
//     int N = sizeof(a) / sizeof(a[0]);
//     int out[N];

//     exe(out, a, b, c, N);

//     for (int i = 0; i < N; i++) {
//         printf("%d ", out[i]);
//     }
//     printf("\n");

//     return 0;
// }

/*


void exe(int *out, int *a, int *b, int *c, int N) {
    svbool_t predicate;
    for (int i = 0; i < N; i+= svcntw()) {
        predicate = svwhilelt_b32(i, N);

        svint32_t vec_a = svld1_s32(predicate, &a[i]);
        svint32_t vec_b = svld1_s32(predicate, &b[i]);
        svint32_t vec_c = svld1_s32(predicate, &c[i]);
        
        svst1_s32(predicate, &out[i], svmul_s32_m(predicate, svadd_s32_m( predicate, vec_a, vec_b), vec_c));
    }
}

*/

/*

void exe(int *out, int *a, int *b, int N) {
    svbool_t predicate;
    for (int i = 0; i < N; i+= svcntw()) {
        predicate = svwhilelt_b32(i, N);

        svint32_t vec_a = svld1_s32(predicate, &a[i]);
        svint32_t vec_b = svld1_s32(predicate, &b[i]);
        
        svst1_s32(predicate, &out[i], svadd_s32_m( predicate, vec_a, vec_b));
    }
}

*/