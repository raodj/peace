#ifndef D2_CUDA_CU
#define D2_CUDA_CU

//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors: Tuan Le                   letm@miamioh.edu
//          Dhananjai M. Rao          raodm@miamioh.edu
//---------------------------------------------------------------------

#include "D2Cuda.h"
#include "ESTCodec.h"
#include "ArgParser.h"
#include "HeuristicChain.h"
#include <algorithm>

// The bitmak to be used when build hash values.
int D2Cuda::BitMask   = 0;
// Instance variable to store the number of bits to be shifted to
// create hash values. This value is initialized to 2*(wordSize-1)
int D2Cuda::bitShift  = 0;

D2Cuda::D2Cuda() : FWAnalyzer("D2Cuda") {
    delta           = NULL;
    alignmentMetric = 0;
    threshold       = 0;
    frameShift      = 1;
}

D2Cuda::~D2Cuda() {
    if (delta != NULL) {
        delete [] delta;
    }
}

void
D2Cuda::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    FWAnalyzer::addCommandLineArguments(argParser);
    // Now add our custom parameters.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--frameShift", "Frame Shift (default=1)",
         &frameShift, ArgParser::INTEGER},
        {"--d2Threshold", "Threshold score to break out of D2Cuda (default=0)",
         &threshold, ArgParser::INTEGER},    
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(ArgsList);
}

bool
D2Cuda::initialize() {
    // Let the base class initialize any additional heuristics
    if (!FWAnalyzer::initialize()) {
        // Error occured when initializing.  This is no good.
        return false;
    }
    // Ensure frameshift is valid.
    if (frameShift < 1) {
        std::cerr << getName()
                  << ": Frame shift must be >= 1"
                  << "(use --frameShift option)\n";
        return false;
    }
    
    // Setup the frequency delta table
    const int MapSize = (1 << (wordSize * 2));
    delta = new int[MapSize];    
    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    BitMask = (1 << (wordSize * 2)) - 1;
    // Compute the number of bits to shift when building hashes
    bitShift = 2 * (wordSize - 1);
    // Set the number of words in a window.
    numWordsInWindow = frameSize - wordSize + 1;    
    // Everything went on well.
    return true;
}

int
D2Cuda::setReferenceEST(const EST* est) {
    ASSERT ( est != NULL );
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(est);
    }
    refEST = est;
    // init ref-est word table
    const char* s1   = est->getSequence();

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << prop.name << std::endl;

    ASSERT("Implement me" == NULL);

    return 1; // For now, return non-zero to signify error
}

float
D2Cuda::getMetric(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    VALIDATE({
        if (otherEST->getID() == refEST->getID()) {
            return 0; // distance to self will be 0
        }
    });
    
    // OK. Run the actual d2 algorithm
    return (float) runD2(otherEST);
}

float
D2Cuda::runD2(const EST* estS2) {
    return 0;  // Currently unimplemented.
}

bool
D2Cuda::getAlignmentData(int &alignmentData) {
    // Simply copy the alignment metric that was computed by the last
    // successful call to the analyze() method.
    alignmentData = alignmentMetric;
    // Let the caller know that the alignment data is available.
    return true;
}

#endif
