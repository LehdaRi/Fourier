#include "Fourier.hpp"

#include <array>
#include <vector>
#include <cmath>
#include <SFML/Graphics.hpp>


using namespace std;
using namespace Fourier;


sf::VertexArray getGraph(const vector<float>& signal,
                         unsigned xstart,
                         unsigned ystart,
                         float xscale,
                         float yscale) {
    const auto s = signal.size();
    sf::VertexArray signalGraph(sf::Lines, (s-1)*2);

    for (auto i=0u; i<s-1; ++i) {
        signalGraph[2*i].position = sf::Vector2f(xstart+xscale*i, ystart-yscale*signal[i]);
        signalGraph[2*i].color = sf::Color(255, 120, 35);
        signalGraph[2*i+1].position = sf::Vector2f(xstart+xscale*(i+1), ystart-yscale*signal[i+1]);
        signalGraph[2*i+1].color = sf::Color(255, 120, 35);
    }

    return signalGraph;
}

sf::VertexArray getGraph(const vector<Compd>& signal,
                         unsigned xstart,
                         unsigned ystart,
                         float xscale,
                         float yscale) {
    const auto s = signal.size();
    sf::VertexArray signalGraph(sf::Lines, (s-1)*4);

    for (auto i=0u; i<s-1; ++i) {
        //  Real part
        signalGraph[2*i].position = sf::Vector2f(xstart+xscale*i, ystart-yscale*signal[i].real());
        signalGraph[2*i].color = sf::Color(255, 120, 35);
        signalGraph[2*i+1].position = sf::Vector2f(xstart+xscale*(i+1), ystart-yscale*signal[i+1].real());
        signalGraph[2*i+1].color = sf::Color(255, 120, 35);
        //  Imaginary part
        signalGraph[2*(s-1)+2*i].position = sf::Vector2f(xstart+xscale*i, ystart-yscale*signal[i].imag());
        signalGraph[2*(s-1)+2*i].color = sf::Color(30, 80, 190);
        signalGraph[2*(s-1)+2*i+1].position = sf::Vector2f(xstart+xscale*(i+1), ystart-yscale*signal[i+1].imag());
        signalGraph[2*(s-1)+2*i+1].color = sf::Color(30, 80, 190);
    }

    return signalGraph;
}


int main(void) {
    const unsigned signalSize = 1024;

    /*sf::RenderWindow window(sf::VideoMode(1024, 800), "SFML window");
    window.setFramerateLimit(60);

    //  Signal
    vector<Compd> signal, transSignal, invSignal;
    signal.resize(signalSize, Compd(0.0f, 0.0f));
    transSignal.resize(signalSize, Compd(0.0f, 0.0f));
    invSignal.resize(signalSize, Compd(0.0f, 0.0f));
    for (auto i=0u; i<signalSize; ++i) {
        signal[i] = Compd((i/64)%2, 0.0f);// powf(cosf(PI2*i*(1.0f/signalSize)) * cosf(PI2*i*(3.0f/signalSize) * sinf(PI2*i*(7.0f/signalSize))), 2.0f));
    }

    Fourier::FFT(&signal[0], &transSignal[0], 1024u, 1u, false, true);
    Fourier::FFT(&transSignal[0], &invSignal[0], 1024u, 1u, true, true);

    auto signalGraph = getGraph(signal, 0, 100, 1.0f, 50.0f);
    auto transSignalGraph = getGraph(transSignal, 0, 300, 1.0f, 0.5f);
    auto invSignalGraph = getGraph(invSignal, 0, 500, 1.0f, 50.0f);

    //  Signal spectra
    vector<float> signalAmpSpectrum;
    vector<float> signalPhaseSpectrum;

    Fourier::AmplitudeSpectrum(transSignal, signalAmpSpectrum, true);
    Fourier::PhaseSpectrum(transSignal, signalPhaseSpectrum, true);

    auto signalAmpSpectrumGraph = getGraph(signalAmpSpectrum, 0, 650, 1.0f, 100.0f);
    auto signalPhaseSpectrumGraph = getGraph(signalPhaseSpectrum, 0, 750, 1.0f, 5.0f);*/

    //  Fourier series
    //float* invSignal = new float[signalSize];
    //Fourier::invDFT(signalAmpSpectrum, signalPhaseSpectrum, invSignal, signalSize, true);
    //auto invSignalGraph = getGraph(invSignal, signalSize, 0, 500, 1.0f, 25.0f);

/*
    {
        printf("Loading image...\n");
        sf::Image srcImg;
        srcImg.loadFromFile("res/perlinNoise2.png");

        const auto srcImgXSize = srcImg.getSize().x;
        const auto srcImgYSize = srcImg.getSize().y;
        const auto* srcImgData = srcImg.getPixelsPtr();

        printf("Converting data...\n");
        std::vector<Compd> emptyData(srcImgXSize*srcImgYSize, Compd(0.0, 0.0));
        std::array<std::vector<Compd>, 4> srcData = { emptyData, emptyData, emptyData, emptyData };
        std::array<std::vector<Compd>, 4> transData = { emptyData, emptyData, emptyData, emptyData };

        for (auto y=0u; y<srcImgYSize; ++y) {
            for (auto x=0u; x<srcImgXSize; ++x) {
                srcData[0][x+y*srcImgXSize] = Compd((double)srcImgData[(x+y*srcImgXSize)*4], 0.0);
                srcData[1][x+y*srcImgXSize] = Compd((double)srcImgData[(x+y*srcImgXSize)*4+1], 0.0);
                srcData[2][x+y*srcImgXSize] = Compd((double)srcImgData[(x+y*srcImgXSize)*4+2], 0.0);
                //srcData[3][x+y*srcImgXSize] = Compd((double)srcImgData[(x+y*srcImgXSize)*4+3], 0.0);
            }
        }

        printf("Transforming...\n");
        printf("Red channel:\n");
        FFT_2D(&srcData[0][0], &transData[0][0], srcImgXSize, srcImgYSize, 1u, false, true);
        printf("Green channel:\n");
        FFT_2D(&srcData[1][0], &transData[1][0], srcImgXSize, srcImgYSize, 1u, false, true);
        printf("Blue channel:\n");
        FFT_2D(&srcData[2][0], &transData[2][0], srcImgXSize, srcImgYSize, 1u, false, true);
        //printf("Alpha channel:\n");
        //FFT_2D(&srcData[3][0], &transData[3][0], srcImgXSize, srcImgYSize, 1u, false, true);

        printf("Calculating amplitude and phase spectra...\n");
        std::vector<double> emptyDataDouble(srcImgXSize*srcImgYSize, 0.0);
        std::array<std::vector<double>, 4> ampSpectrum { emptyDataDouble, emptyDataDouble, emptyDataDouble, emptyDataDouble };
        std::array<std::vector<double>, 4> phaseSpectrum { emptyDataDouble, emptyDataDouble, emptyDataDouble, emptyDataDouble };

        AmplitudeSpectrum(transData[0], ampSpectrum[0], true);
        AmplitudeSpectrum(transData[1], ampSpectrum[1], true);
        AmplitudeSpectrum(transData[2], ampSpectrum[2], true);
        //AmplitudeSpectrum(transData[3], ampSpectrum[3], true);
        PhaseSpectrum(transData[0], phaseSpectrum[0], true);
        PhaseSpectrum(transData[1], phaseSpectrum[1], true);
        PhaseSpectrum(transData[2], phaseSpectrum[2], true);
        //PhaseSpectrum(transData[3], phaseSpectrum[3], true);

        printf("Writing amplitude and phase spectra in files...\n");
        sf::Image ampSpectrumImg, phaseSpectrumImg;
        ampSpectrumImg.create(srcImgXSize, srcImgYSize);
        phaseSpectrumImg.create(srcImgXSize, srcImgYSize);

        for (auto y=0u; y<srcImgYSize; ++y) {
            for (auto x=0u; x<srcImgXSize; ++x) {
                auto xx = (x+srcImgXSize/2)%srcImgXSize;
                auto yy = (y+srcImgYSize/2)%srcImgYSize;

                double r = 255*(0.1*std::log(ampSpectrum[0][x+y*srcImgXSize]/255.0)+1.0);
                if (r<0.0) r = 0.0;
                double g = 255*(0.1*std::log(ampSpectrum[1][x+y*srcImgXSize]/255.0)+1.0);
                if (g<0.0) g = 0.0;
                double b = 255*(0.1*std::log(ampSpectrum[2][x+y*srcImgXSize]/255.0)+1.0);
                if (b<0.0) b = 0.0;
                double a = 255.0;
                //double a = 255*(0.1*std::log(ampSpectrum[3][x+y*srcImgXSize]/255.0)+1.0);
                //if (a<0.0) a = 0.0;

                ampSpectrumImg.setPixel(xx, yy, sf::Color(r, g, b, a));

                phaseSpectrumImg.setPixel(xx, yy, sf::Color(255*(phaseSpectrum[0][x+y*srcImgXSize]/PI2),
                                                            255*(phaseSpectrum[1][x+y*srcImgXSize]/PI2),
                                                            255*(phaseSpectrum[2][x+y*srcImgXSize]/PI2)));//,
                                                            //255*(phaseSpectrum[3][x+y*srcImgXSize]/PI2)));
            }
        }

        ampSpectrumImg.saveToFile("ampSpectrum.png");
        phaseSpectrumImg.saveToFile("phaseSpectrum.png");
    }
*/

    {
        printf("Inverse transforming from spectra...\n");
        sf::Image ampSpectrumImg, phaseSpectrumImg;
        ampSpectrumImg.loadFromFile("ampSpectrum.png");
        phaseSpectrumImg.loadFromFile("phaseSpectrum.png");

        const auto spectraImgXSize = ampSpectrumImg.getSize().x;
        const auto spectraImgYSize = ampSpectrumImg.getSize().y;
        const auto* ampSpectrumImgData = ampSpectrumImg.getPixelsPtr();
        const auto* phaseSpectrumImgData = phaseSpectrumImg.getPixelsPtr();

        printf("Converting data...\n");
        std::vector<Compd> emptyDataCompd(spectraImgXSize*spectraImgYSize, Compd(0.0, 0.0));
        std::array<std::vector<Compd>, 4> spectrumData { emptyDataCompd, emptyDataCompd, emptyDataCompd, emptyDataCompd };
        for (auto y=0u; y<spectraImgYSize; ++y) {
            for (auto x=0u; x<spectraImgXSize; ++x) {
                auto xx = (x+spectraImgXSize/2)%spectraImgXSize;
                auto yy = (y+spectraImgYSize/2)%spectraImgYSize;

                //for (auto i=0u; i<4; ++i) {
                for (auto i=0u; i<3; ++i) {
                    double amp = 255.0*std::exp((10.0*(ampSpectrumImgData[(xx+yy*spectraImgXSize)*4+i]/255.0)-10.0));
                    amp *= spectraImgXSize*spectraImgYSize;
                    double phase = phaseSpectrumImgData[(xx+yy*spectraImgXSize)*4+i]*(PI2/255.0);
                    spectrumData[i][x+y*spectraImgXSize] =  Compd(amp*std::cos(phase), amp*std::sin(phase));
                }
            }
        }

        printf("Transforming...\n");
        std::array<std::vector<Compd>, 4> spectraInvData { emptyDataCompd, emptyDataCompd, emptyDataCompd, emptyDataCompd };
        FFT_2D(&spectrumData[0][0], &spectraInvData[0][0], spectraImgXSize, spectraImgYSize, 1u, true, true);
        FFT_2D(&spectrumData[1][0], &spectraInvData[1][0], spectraImgXSize, spectraImgYSize, 1u, true, true);
        FFT_2D(&spectrumData[2][0], &spectraInvData[2][0], spectraImgXSize, spectraImgYSize, 1u, true, true);
        //FFT_2D(&spectrumData[3][0], &spectraInvData[3][0], spectraImgXSize, spectraImgYSize, 1u, true, true);

        sf::Image invImg;
        invImg.create(spectraImgXSize, spectraImgYSize);
        for (auto y=0u; y<spectraImgYSize; ++y) {
            for (auto x=0u; x<spectraImgXSize; ++x) {
                int r = spectraInvData[0][x+y*spectraImgXSize].real();
                if (r<0) r=0; else if (r>255) r=255;
                int g = spectraInvData[1][x+y*spectraImgXSize].real();
                if (g<0) g=0; else if (g>255) g=255;
                int b = spectraInvData[2][x+y*spectraImgXSize].real();
                if (b<0) b=0; else if (b>255) b=255;
                //int a = spectraInvData[3][x+y*spectraImgXSize].real();
                //if (a<0) a=0; else if (a>255) a=255;
                int a = 255;

                invImg.setPixel(x, y, sf::Color(r, g, b, a));
            }
        }
        invImg.saveToFile("spectraInverse.png");
    }

/*
    while (window.isOpen())
    {
        // Event processing
        sf::Event event;
        while (window.pollEvent(event))
        {
           // Request for closing the window
           if (event.type == sf::Event::Closed)
               window.close();
        }

        window.clear();

        window.draw(signalGraph);
        window.draw(transSignalGraph);
        window.draw(invSignalGraph);
        window.draw(signalAmpSpectrumGraph);
        window.draw(signalPhaseSpectrumGraph);
        window.display();
    }*/
}
