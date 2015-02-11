#ifndef FOURIER_HPP
#define FOURIER_HPP


#include <vector>
#include <complex>


#define PI2 6.28318530718
#define I() Compd(0.0, 1.0)


namespace Fourier {

    typedef std::complex<float> Compf;
    typedef std::complex<double> Compd;

    void DFT(const std::vector<Compd>& src, std::vector<Compd>& dest, bool inverse = false, bool printTime = false);

    void AmplitudeSpectrum(const std::vector<Compd>& src, std::vector<double>& dest, bool printTime = false);
    void PhaseSpectrum(const std::vector<Compd>& src, std::vector<double>& dest, bool printTime = false);

    void FFT(const Compd* src, Compd* dest, unsigned size, unsigned stride = 1u, bool inverse = false, bool printTime = false);
    void FFT_recursion(const Compd* src, Compd* dest, unsigned size, unsigned stride, bool inverse);

    void FFT_2D(const Compd* src, Compd* dest, unsigned xsize, unsigned ysize, unsigned stride = 1u, bool inverse = false, bool printTime = false);

    void series(const double* amp, const double* phase, double* dest, unsigned size, bool printTime);

} // namespace Fourier


#endif // FOURIER_HPP
