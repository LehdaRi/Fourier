#include "Fourier.hpp"

#include <cmath>
#include <chrono>
#include <iostream>


void Fourier::DFT(const float* src, float* amp, float* phase, unsigned size, bool printTime) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for (auto k=0u; k<size; ++k) {
        float re(0.0f), img(0.0f);

        for (auto n=0u; n<size; ++n) {
            float a = (float)(PI2*k*n)/size;
            re += src[n]*cosf(a);
            img += src[n]*sinf(-a);
        }

        amp[k] = sqrtf(powf(re, 2.0f) + powf(img, 2.0f)) / size;
        phase[k] = atan2f(img, re);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Calculating DFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
}

void Fourier::FFT(const float* src, float* amp, float* phase, unsigned size, bool printTime, unsigned stride) {
    std::chrono::steady_clock::time_point start, end;

    if (printTime)
        start = std::chrono::steady_clock::now();



    if (size == 1) {
        amp[0] = src[0];
        phase[0] = 0.0f;
        return;
    }

    FFT(src, amp, phase, size/2, false, 2*stride);
    FFT(src+stride, amp+size/2, phase+size/2, size/2, false, 2*stride);

    float re(0.0f), img(0.0f);

    for (auto k=0u; k<size/2; ++k) {
        float a = (float)(PI2*k)/size;
        float t = src[k];

        re = t + src[k+size/2]*cosf(a);
        img = t + src[k+size/2]*sinf(-a);
        amp[k] = sqrtf(powf(re, 2.0f) + powf(img, 2.0f)) / size;
        phase[k] = atan2f(img, re);

        re = t - src[k+size/2]*cosf(a);
        img = t - src[k+size/2]*sinf(-a);
        amp[k+size/2] = sqrtf(powf(re, 2.0f) + powf(img, 2.0f)) / size;
        phase[k+size/2] = atan2f(img, re);
    }

    if (printTime) {
        end = std::chrono::steady_clock::now();
        std::cout << "Calculating FFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    }
}

void Fourier::invDFT(const float* amp, const float* phase, float* dest, unsigned size, bool printTime = false) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for(auto k = 0u; k < size; k++) {
        for(auto n = 0u; n < size; n++) {
            dest[k] += amp[n]*cosf(PI2*k*(float)n/size+phase[n]);
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Calculating fourier series took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
}

/*void Fourier::dft2d(const float* src, float* amp, float* phase, unsigned size) {
    for (auto k=0u; k<size; ++k) {
        float re(0.0f), img(0.0f);

        for (auto n=0u; n<size; ++n) {
            float a = (float)(PI2*k*n)/size;
            re += src[n]*cosf(a);
            img += src[n]*sinf(-a);
        }

        amp[k] = sqrtf(powf(re, 2.0f) + powf(img, 2.0f)) / size;
        phase[k] = atan2f(img, re);
    }
}*/
