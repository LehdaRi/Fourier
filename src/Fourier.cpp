#include "Fourier.hpp"

#include <cmath>
#include <chrono>
#include <iostream>


void Fourier::DFT(const std::vector<Compd>& src, std::vector<Compd>& dest, bool inverse, bool printTime) {
    std::chrono::steady_clock::time_point start, end;
    if (printTime)
        start = std::chrono::steady_clock::now();

    dest.resize(src.size(), Compd(0.0, 0.0));

    double a;
    if (inverse)
        a = PI2/src.size();
    else
        a = -PI2/src.size();

    for (auto k=0u; k<src.size(); ++k) {
        for (auto n=0u; n<src.size(); ++n) {
            dest[k] += src[n]*std::exp(Compd(0.0, 1.0)*(a*k*n));
        }

        if (inverse)
            dest[k] /= (double)src.size();
    }

    if (printTime) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Calculating DFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    }
}

void Fourier::AmplitudeSpectrum(const std::vector<Compd>& src, std::vector<double>& dest, bool printTime) {
    std::chrono::steady_clock::time_point start, end;
    if (printTime)
        start = std::chrono::steady_clock::now();

    dest.clear();
    //dest.reserve(src.size());

    /*for (auto k=0u; k<size; ++k) {
        if (dest[k].real() > 0.0 || dest[k].imag() > 0.0)
            printf("dest[%u] = %0.2f %c %0.2fi\n", k, dest[k].real(), dest[k].imag() < 0 ? '-' : '+', dest[k].imag());
    }*/

    for (auto& s : src) {
        dest.push_back(std::sqrt(std::pow(s.real(), 2.0) + std::pow(s.imag(), 2.0)) / src.size());
        //if (s.real() > 0.0 || s.imag() > 0.0)
            //printf("%0.2f%s%0.2fi\n", s.real(), s.imag() >= 0 ? " +" : " ", s.imag());
    }

    if (printTime) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Calculating amplitude spectrum took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    }
}

void Fourier::PhaseSpectrum(const std::vector<Compd>& src, std::vector<double>& dest, bool printTime) {
    std::chrono::steady_clock::time_point start, end;
    if (printTime)
        start = std::chrono::steady_clock::now();

    dest.clear();
    dest.reserve(src.size());

    for (auto& s : src)
        dest.push_back(std::atan2(s.imag(), s.real()));

    if (printTime) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Calculating phase spectrum took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    }
}

void Fourier::FFT(const Compd* src, Compd* dest, unsigned size, unsigned stride, bool inverse, bool printTime) {
    std::chrono::steady_clock::time_point start, end;
    if (printTime)
        start = std::chrono::steady_clock::now();

    const unsigned hs = size/2;

    double a;
    if (inverse)
        a = PI2/size;
    else
        a = -PI2/size;

    if (size == 1) {
        dest[0] = src[0];
    }
    else {
        FFT_recursion(src, dest, hs, 2*stride, inverse);
        FFT_recursion(src+stride, dest+hs, hs, 2*stride, inverse);

        for (auto k=0u; k<hs; ++k) {
            auto temp = dest[k];
            dest[k] = temp + dest[k+hs]*std::exp(I()*(a*k));
            dest[k+hs] = temp - dest[k+hs]*std::exp(I()*(a*k));
        }

        if (inverse) {
            for (auto k=0u; k<size; ++k)
                dest[k] /= (double)size;
        }
    }

    if (printTime) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Calculating FFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    }
}

void Fourier::FFT_recursion(const Compd* src, Compd* dest, unsigned size, unsigned stride, bool inverse) {
    const unsigned hs = size/2;

    double a;
    if (inverse)
        a = PI2/size;
    else
        a = -PI2/size;

    if (size == 1) {
        dest[0] = src[0];
    }
    else {
        FFT_recursion(src, dest, hs, 2*stride, inverse);
        FFT_recursion(src+stride, dest+hs, hs, 2*stride, inverse);

        for (auto k=0u; k<hs; ++k) {
            auto temp = dest[k];
            dest[k] = temp + dest[k+hs]*std::exp(I()*(a*k));
            dest[k+hs] = temp - dest[k+hs]*std::exp(I()*(a*k));
        }
    }
}

void Fourier::FFT_2D(const Compd* src, Compd* dest, unsigned xsize, unsigned ysize, unsigned stride, bool inverse, bool printTime) {
    std::chrono::steady_clock::time_point start, end;
    if (printTime)
        start = std::chrono::steady_clock::now();

    std::vector<Compd> tempSrc(ysize, Compd(0.0, 0.0));
    std::vector<Compd> tempDest(ysize, Compd(0.0, 0.0));

    for (auto y=0u; y<ysize; ++y) {
        FFT(src+xsize*y, &dest[xsize*y], xsize, 1u, inverse, false);
    }

    for (auto x=0u; x<xsize; ++x) {
        for (auto y=0u; y<ysize; ++y)
            tempSrc[y] = dest[x+y*xsize];

        FFT(&tempSrc[0], &tempDest[0], ysize, 1u, inverse, false);

        for (auto y=0u; y<ysize; ++y)
            dest[x+y*xsize] = tempDest[y];
    }

    /*for (auto y=0u; y<ysize; ++y)
        for (auto x=0u; x<xsize; ++x)
            dest[x+y*xsize] = temp[x+y*xsize];*/

    if (printTime) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Calculating 2D FFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    }
}

void Fourier::series(const double* amp, const double* phase, double* dest, unsigned size, bool printTime = false) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for(auto k = 0u; k < size; k++) {
        for(auto n = 0u; n < size; n++) {
            dest[k] += amp[n]*cos(PI2*k*(double)n/size+phase[n]);
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Calculating fourier series took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
}
