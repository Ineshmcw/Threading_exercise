#include <iostream>
#include <vector>
#include <arm_sve.h>

// Convert an RGB image to grayscale using SVE intrinsics
void rgb_to_grayscale_sve(uint8_t* rgb, uint8_t* gray, size_t num_pixels) {
    size_t i = 0;
    svfloat32_t coeff_r = svdup_f32(0.2989f);
    svfloat32_t coeff_g = svdup_f32(0.5870f);
    svfloat32_t coeff_b = svdup_f32(0.1140f);

    size_t vec_size = svcntb() / 3; // Number of RGB pixels per iteration

    for (; i + vec_size <= num_pixels; i += vec_size) {
        svbool_t pg = svwhilelt_b8(i, num_pixels);
        svuint8x3_t rgb_vec = svld3_u8(pg, &rgb[i * 3]);

        svuint8_t r = svget3(rgb_vec, 0);
        svuint8_t g = svget3(rgb_vec, 1);
        svuint8_t b = svget3(rgb_vec, 2);

        // Convert uint8 to float32
        svfloat32_t r_f32 = svcvt_f32_u32_z(pg, svreinterpret_u32_u8(r));
        svfloat32_t g_f32 = svcvt_f32_u32_z(pg, svreinterpret_u32_u8(g));
        svfloat32_t b_f32 = svcvt_f32_u32_z(pg, svreinterpret_u32_u8(b));

        // Apply grayscale formula
        svfloat32_t gray_f32 = svmla_f32_z(pg, svmla_f32_z(pg, svmul_f32_z(pg, r_f32, coeff_r), g_f32, coeff_g), b_f32, coeff_b);

        // Convert float32 back to uint32
        svuint32_t gray_u32 = svcvt_u32_f32_z(pg, gray_f32);

        // Narrow from uint32 -> uint16 -> uint8
        svuint16_t gray_u16 = svqxtnt_u32(svundef_u16(), gray_u32);
        svuint8_t gray_u8 = svqxtnt_u16(svundef_u8(), gray_u16);

        // Store the grayscale values
        svst1_u8(pg, &gray[i], gray_u8);
    }

    // Process remaining pixels
    for (; i < num_pixels; i++) {
        uint8_t r = rgb[i * 3];
        uint8_t g = rgb[i * 3 + 1];
        uint8_t b = rgb[i * 3 + 2];
        gray[i] = static_cast<uint8_t>((r * 0.2989f) + (g * 0.5870f) + (b * 0.1140f));
    }
}

int main() {
    size_t num_pixels = 1024; // Example image size
    std::vector<uint8_t> rgb(num_pixels * 3, 128); // Dummy RGB data
    std::vector<uint8_t> gray(num_pixels, 0);

    rgb_to_grayscale_sve(rgb.data(), gray.data(), num_pixels);

    std::cout << "Grayscale conversion using SVE completed.\n";
    return 0;
}
