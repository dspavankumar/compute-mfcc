// -----------------------------------------------------------------------------
// Wrapper for MFCC feature extractor
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

#include <iostream>
#include <fstream>
#include "mfcc.cpp"

int main(int argc, char* argv[]) {
    // Check arguments
    if (argc != 3) {
        std::cout << "Usage: compute-mfcc input.wav output.mfc" << std::endl;
        return 0;
    }
    
    // Assign paths
    const char *wavPath=argv[1], *mfcPath=argv[2];

    // Initialise input and output streams    
    std::ifstream wavFp;
    wavFp.open(wavPath);
    std::ofstream mfcFp;
    mfcFp.open(mfcPath);
    
    // Check if input is readable
    if (!wavFp.is_open()) {
        std::cerr << "Unable to open input file: " << wavPath << std::endl;
        return 1;
    }
    // Check if output is writable
    if (!mfcFp.is_open()) {
        std::cerr << "Unable to open output file: " << mfcPath << std::endl;
        wavFp.close();
        return 1;
    }

    // Initialise MFCC class instance
    int winLength=25, frameShift=10; // Both in milliseconds
    int fs=16000;        // Sampling frequency in Hertz
    int numCepstra = 12; // We get log-energy extra.
    MFCC mfccComputer (fs, numCepstra, winLength, frameShift);
    
    // Process
    if (mfccComputer.process (wavFp, mfcFp))
        std::cerr << "Error processing " << wavPath << std::endl;

    wavFp.close();
    mfcFp.close();
    return 0;
}
