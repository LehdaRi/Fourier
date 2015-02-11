#include "Fourier.hpp"

#include <vector>
#include <cmath>
#include <SFML/Graphics.hpp>


#define E 2.718281828459045


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

    //  Signal spectri
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



    sf::Image srcImg;
    printf("Loading image...\n");
    srcImg.loadFromFile("res/fractal.png");

    const auto srcImgXSize = srcImg.getSize().x;
    const auto srcImgYSize = srcImg.getSize().y;
    const auto* srcImgData = srcImg.getPixelsPtr();

    printf("Converting data...\n");
    std::vector<Compd> srcData(srcImgXSize*srcImgYSize, Compd(0.0, 0.0));
    std::vector<Compd> transData(srcImgXSize*srcImgYSize, Compd(0.0, 0.0));

    for (auto y=0u; y<srcImgYSize; ++y)
        for (auto x=0u; x<srcImgXSize; ++x) {
            srcData[x+y*srcImgXSize] = Compd((double)srcImgData[(x+y*srcImgXSize)*4], 0.0);
        }

    printf("Transforming...\n");
    FFT_2D(&srcData[0], &transData[0], srcImgXSize, srcImgYSize, 1u, false, true);

    for (auto i=0u; i<10; ++i)
        printf("transData[%u]: %0.2f %0.2fi\n", i, transData[i].real(), transData[i].imag());

    /*{
        printf("Inverse transforming...\n");
        std::vector<Compd> invData(srcImgXSize*srcImgYSize, Compd(0.0, 0.0));
        FFT_2D(&transData[0], &invData[0], srcImgXSize, srcImgYSize, 1u, true, true);
        sf::Image invImg;
        invImg.create(srcImgXSize, srcImgYSize);
        for (auto y=0u; y<srcImgYSize; ++y) {
            for (auto x=0u; x<srcImgXSize; ++x) {
                invImg.setPixel(x, y, sf::Color(invData[x+y*srcImgXSize].real(), 0, 0));
            }
        }
        invImg.saveToFile("inverse.png");
    }*/

    std::vector<double> ampSpectrum(srcImgXSize*srcImgYSize, 0.0);
    std::vector<double> phaseSpectrum(srcImgXSize*srcImgYSize, 0.0);

    AmplitudeSpectrum(transData, ampSpectrum, true);
    PhaseSpectrum(transData, phaseSpectrum, true);

    printf("Writing amplitude and phase spectri in files...\n");
    sf::Image ampSpectrumImg, phaseSpectrumImg;
    ampSpectrumImg.create(srcImgXSize, srcImgYSize);
    phaseSpectrumImg.create(srcImgXSize, srcImgYSize);

    for (auto i=0u; i<10; ++i)
        printf("ampSpectrum[%u]: %0.2f\n", i, ampSpectrum[i]);


    for (auto y=0u; y<srcImgYSize; ++y) {
        for (auto x=0u; x<srcImgXSize; ++x) {
            auto xx = (x+srcImgXSize/2)%srcImgXSize;
            auto yy = (y+srcImgYSize/2)%srcImgYSize;
            double r = 255*(0.1*std::log(ampSpectrum[x+y*srcImgXSize]/255.0)+1.0);
            if (r<0.0) r = 0.0;
            ampSpectrumImg.setPixel(xx, yy, sf::Color(r, 0, 0));
            phaseSpectrumImg.setPixel(xx, yy, sf::Color(255*(phaseSpectrum[x+y*srcImgXSize]/PI2), 0, 0));
        }
    }

    ampSpectrumImg.saveToFile("ampSpectrum.png");
    phaseSpectrumImg.saveToFile("phaseSpectrum.png");

    printf("Inverse transforming from spectri...\n");
    {
        //sf::Image ampSpectrumImg, phaseSpectrumImg;
        //ampSpectrumImg.loadFromFile(fileN)

        const auto spectriImgXSize = ampSpectrumImg.getSize().x;
        const auto spectriImgYSize = ampSpectrumImg.getSize().y;
        const auto* ampSpectrumImgData = ampSpectrumImg.getPixelsPtr();
        const auto* phaseSpectrumImgData = phaseSpectrumImg.getPixelsPtr();

        std::vector<Compd> spectrumData(spectriImgXSize*spectriImgYSize, Compd(0.0, 0.0));
        for (auto y=0u; y<spectriImgYSize; ++y) {
            for (auto x=0u; x<spectriImgXSize; ++x) {
                auto xx = (x+srcImgXSize/2)%srcImgXSize;
                auto yy = (y+srcImgYSize/2)%srcImgYSize;
                double amp = 255.0*std::exp((10.0*(ampSpectrumImgData[(xx+yy*spectriImgXSize)*4]/255.0)-10.0));
                amp *= spectriImgXSize*spectriImgYSize;
                double phase = phaseSpectrumImgData[(xx+yy*spectriImgXSize)*4]*(PI2/255.0);
                spectrumData[x+y*spectriImgXSize] =  Compd(amp*std::cos(phase), amp*std::sin(phase));
            }
        }

        for (auto i=0u; i<10; ++i)
            printf("spectrumData[%u]: %0.2f %0.2fi\n", i, spectrumData[i].real(), spectrumData[i].imag());

        std::vector<Compd> spectriInvData(spectriImgXSize*spectriImgYSize, Compd(0.0, 0.0));
        FFT_2D(&spectrumData[0], &spectriInvData[0], spectriImgXSize, spectriImgYSize, 1u, true, true);

        sf::Image invImg;
        invImg.create(spectriImgXSize, spectriImgYSize);
        for (auto y=0u; y<spectriImgYSize; ++y) {
            for (auto x=0u; x<spectriImgXSize; ++x) {
                int r = spectriInvData[x+y*spectriImgXSize].real();
                if (r<0) r=0;
                if (r>255) r=255;
                invImg.setPixel(x, y, sf::Color(r, 0, 0));
            }
        }
        invImg.saveToFile("spectriInverse.png");
    }

    /*while (window.isOpen())
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
