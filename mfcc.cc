// -----------------------------------------------------------------------------
//  A simple MFCC extractor using C++ STL and C++11
// -----------------------------------------------------------------------------
//
//  Copyright (C) 2016 D S Pavan Kumar
//  dspavankumar [at] gmail [dot] com
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include<algorithm>
#include<numeric>
#include<complex>
#include<vector>
#include<map>
#include<math.h>
#include"wavHeader.h"

typedef std::vector<double> v_d;
typedef std::complex<double> c_d;
typedef std::vector<v_d> m_d;
typedef std::vector<c_d> v_c_d;
typedef std::map<int,std::map<int,c_d> > twmap;

class MFCC {

private:
    const double PI = 4*atan(1.0);   // Pi = 3.14...
    int fs;
    twmap twiddle;
    size_t winLengthSamples, frameShiftSamples, numCepstra, numFFT, numFFTBins, numFilters;
    double preEmphCoef, lowFreq, highFreq;
    v_d frame, powerSpectralCoef, lmfbCoef, hamming, mfcc, prevsamples;
    m_d fbank, dct;

private:
    // Hertz to Mel conversion
    inline double hz2mel (double f) {
        return 2595*std::log10 (1+f/700);
    }

    // Mel to Hertz conversion
    inline double mel2hz (double m) {
        return 700*(std::pow(10,m/2595)-1);
    }

    // Twiddle factor computation
    void compTwiddle(void) {
        const c_d J(0,1);      // Imaginary number 'j'
        for (int N=2; N<=numFFT; N*=2)
            for (int k=0; k<=N/2-1; k++)
                twiddle[N][k] = exp(-2*PI*k/N*J);
    }
    
    // Cooley-Tukey DIT-FFT recursive function
    v_c_d fft(v_c_d x) {
        int N = x.size();
        if (N==1)
            return x;
    
        v_c_d xe(N/2,0), xo(N/2,0), Xjo, Xjo2;
        int i;
    
        // Construct arrays from even and odd indices
        for (i=0; i<N; i+=2)
            xe[i/2] = x[i];
        for (i=1; i<N; i+=2)
            xo[(i-1)/2] = x[i];
    
        // Compute N/2-point FFT
        Xjo = fft(xe);
        Xjo2 = fft(xo);
        Xjo.insert (Xjo.end(), Xjo2.begin(), Xjo2.end());
    
        // Butterfly computations
        for (i=0; i<=N/2-1; i++) {
            c_d t = Xjo[i], tw = twiddle[N][i];
            Xjo[i] = t + tw * Xjo[i+N/2];
            Xjo[i+N/2] = t - tw * Xjo[i+N/2];
        }
        return Xjo;
    }

    //// Frame processing routines
    // Pre-emphasis and Hamming window
    void preEmphHam(void) {
        v_d procFrame(frame.size(), hamming[0]*frame[0]);
        for (int i=1; i<frame.size(); i++)
            procFrame[i] = hamming[i] * (frame[i] - preEmphCoef * frame[i-1]);
        frame = procFrame;
    }

    // Power spectrum computation
    void computePowerSpec(void) {
        frame.resize(numFFT); // Pads zeros
        v_c_d framec (frame.begin(), frame.end()); // Complex frame
        v_c_d fftc = fft(framec);
    
        for (int i=0; i<numFFTBins; i++)
            powerSpectralCoef[i] = pow(abs(fftc[i]),2);
    }

    // Applying log Mel filterbank (LMFB)
    void applyLMFB(void) {
        lmfbCoef.assign(numFilters,0);
        
        for (int i=0; i<numFilters; i++) {
            // Multiply the filterbank matrix
            for (int j=0; j<fbank[i].size(); j++)
                lmfbCoef[i] += fbank[i][j] * powerSpectralCoef[j];
            // Apply Mel-flooring
            if (lmfbCoef[i] < 1.0)
                lmfbCoef[i] = 1.0;
        }
        
        // Applying log on amplitude
        for (int i=0; i<numFilters; i++)
            lmfbCoef[i] = std::log (lmfbCoef[i]);
    }
    
    // Computing discrete cosine transform
    void applyDct(void) {
        mfcc.assign(numCepstra+1,0);
        for (int i=0; i<=numCepstra; i++) {
            for (int j=0; j<numFilters; j++)
                mfcc[i] += dct[i][j] * lmfbCoef[j];
        }
    }
   
    // Initialisation routines
    // Pre-computing Hamming window and dct matrix
    void initHamDct(void) {
        int i, j;

        hamming.assign(winLengthSamples,0);
        for (i=0; i<winLengthSamples; i++)
            hamming[i] = 0.54 - 0.46 * cos(2 * PI * i / (winLengthSamples-1));

        v_d v1(numCepstra+1,0), v2(numFilters,0);
        for (i=0; i <= numCepstra; i++)
            v1[i] = i;
        for (i=0; i < numFilters; i++)
            v2[i] = i + 0.5;

        dct.reserve (numFilters*(numCepstra+1));        
        double c = sqrt(2.0/numFilters);
        for (i=0; i<=numCepstra; i++) {
            v_d dtemp;
            for (j=0; j<numFilters; j++)
                dtemp.push_back (c * cos(PI / numFilters * v1[i] * v2[j]));
            dct.push_back(dtemp);
        }
    }

    // Precompute filterbank
    void initFilterbank () {
        // Convert low and high frequencies to Mel scale
        double lowFreqMel = hz2mel(lowFreq);
        double highFreqMel = hz2mel (highFreq);

        // Calculate filter centre-frequencies
        v_d filterCentreFreq;
        filterCentreFreq.reserve (numFilters+2);
        for (int i=0; i<numFilters+2; i++)
            filterCentreFreq.push_back (mel2hz(lowFreqMel + (highFreqMel-lowFreqMel)/(numFilters+1)*i));

        // Calculate FFT bin frequencies
        v_d fftBinFreq;
        fftBinFreq.reserve(numFFTBins);
        for (int i=0; i<numFFTBins; i++)
            fftBinFreq.push_back (fs/2.0/(numFFTBins-1)*i);
            
        // Filterbank: Allocate memory
        fbank.reserve (numFilters*numFFTBins);
        
        // Populate the fbank matrix
        for (int filt=1; filt<=numFilters; filt++) {
            v_d ftemp;
            for (int bin=0; bin<numFFTBins; bin++) {
                double weight;
                if (fftBinFreq[bin] < filterCentreFreq[filt-1])
                    weight = 0;
                else if (fftBinFreq[bin] <= filterCentreFreq[filt])
                    weight = (fftBinFreq[bin] - filterCentreFreq[filt-1]) / (filterCentreFreq[filt] - filterCentreFreq[filt-1]);
                else if (fftBinFreq[bin] <= filterCentreFreq[filt+1])
                    weight = (filterCentreFreq[filt+1] - fftBinFreq[bin]) / (filterCentreFreq[filt+1] - filterCentreFreq[filt]);
                else
                    weight = 0;
                ftemp.push_back (weight);
            }
            fbank.push_back(ftemp);
        }
    }

    // Convert vector of double to string (for writing MFCC file output)
    std::string v_d_to_string (v_d vec) {
        std::stringstream vecStream;
        for (int i=0; i<vec.size()-1; i++) {
            vecStream << std::scientific << vec[i];
            vecStream << ", ";
        }
        vecStream << std::scientific << vec.back();
        vecStream << "\n";
        return vecStream.str();
    }

public:
    // MFCC class constructor
    MFCC(int sampFreq=16000, int nCep=12, int winLength=25, int frameShift=10, int numFilt=40, double lf=50, double hf=6500) {
        fs          = sampFreq;     // Sampling frequency
        numCepstra  = nCep;         // Number of cepstra
        numFFT      = 512;          // FFT size
        numFilters  = numFilt;      // Number of Mel warped filters
        preEmphCoef = 0.97;         // Pre-emphasis coefficient
        lowFreq     = lf;           // Filterbank low frequency cutoff in Hertz
        highFreq    = hf;         // Filterbank high frequency cutoff in Hertz
        winLengthSamples   = winLength * fs / 1e3;  // winLength in milliseconds
        frameShiftSamples  = frameShift * fs / 1e3; // frameShift in milliseconds
        
        numFFTBins = numFFT/2 + 1;
        powerSpectralCoef.assign (numFFTBins, 0);
        prevsamples.assign (winLengthSamples-frameShiftSamples, 0);

        initFilterbank();
        initHamDct();
        compTwiddle();
    }

    // Process each frame and extract MFCC
    std::string processFrame(int16_t* samples, size_t N) {
        // Add samples from the previous frame that overlap with the current frame
        // to the current samples and create the frame.
        frame = prevsamples;
        for (int i=0; i<N; i++)
            frame.push_back(samples[i]);
        prevsamples.assign(frame.begin()+frameShiftSamples, frame.end());

        preEmphHam();
        computePowerSpec();
        applyLMFB();
        applyDct();

        return v_d_to_string (mfcc);
    }

    // Read input file stream, extract MFCCs and write to output file stream
    int process (std::ifstream &wavFp, std::ofstream &mfcFp) {
        // Read the wav header    
        wavHeader hdr;
        int headerSize = sizeof(wavHeader);
        wavFp.read((char *) &hdr, headerSize);

        // Check audio format
        if (hdr.AudioFormat != 1 || hdr.bitsPerSample != 16) {
            std::cerr << "Unsupported audio format, use 16 bit PCM Wave" << std::endl;
            return 1;
        }
        // Check sampling rate
        if (hdr.SamplesPerSec != fs) {
            std::cerr << "Sampling rate mismatch: Found " << hdr.SamplesPerSec << " instead of " << fs <<std::endl;
            return 1;
        }
        
        // Initialise buffer
        uint16_t bufferLength = winLengthSamples-frameShiftSamples;
        int16_t* buffer = new int16_t[bufferLength];
        int bufferBPS = (sizeof buffer[0]);

        // Read and set the initial samples        
        wavFp.read((char *) buffer, bufferLength*bufferBPS);
        for (int i=0; i<bufferLength; i++)
            prevsamples[i] = buffer[i];        
        delete [] buffer;
        
        // Recalculate buffer size
        bufferLength = frameShiftSamples;
        buffer = new int16_t[bufferLength];
        
        // Read data and process each frame
        wavFp.read((char *) buffer, bufferLength*bufferBPS);
        while (wavFp.gcount() == bufferLength*bufferBPS && !wavFp.eof()) {
            mfcFp << processFrame(buffer, bufferLength);
            wavFp.read((char *) buffer, bufferLength*bufferBPS);
        }
        delete [] buffer;
        buffer = nullptr;
        return 0;
    }
};
