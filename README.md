## A Simple MFCC Feature Extractor using C++ STL and C++11
# Features
* Takes PCM Wave input.
* Outputs MFCCs as comma separated floating point values, each line representing a frame.
* Applies pre-emphasis, Hamming window, computes FFT and power spectrum.
* Applies equal height log-Mel filterbank followed by discrete cosine transform.
* Window length, frame shift, sampling rate and number of cepstra are configurable in the wrapper script.
* Outputs log-energy as the first coefficient, followed by the specified number of cepstra.

# Compilation
g++ -std=c++11 -O3 compute-mfcc.cc -o compute-mfcc

# Usage
./compute-mfcc input.wav output.mfc

# License
GNU GPL V3.0

# Contributors
D S Pavan Kumar
