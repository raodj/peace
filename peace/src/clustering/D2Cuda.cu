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

D2Cuda::D2Cuda() : FWAnalyzer("D2Cuda") {
    alignmentMetric  = 0;
    threshold        = 0;
    frameShift       = 1;
    // Fix frame and word size
    frameSize        = 100;
    wordSize         = 6;
    bitMask          = (1 << (wordSize * 2)) - 1; // 4095
    numWordsInWindow = frameSize - wordSize + 1;
}

D2Cuda::~D2Cuda() {

}

void D2Cuda::addCommandLineArguments(ArgParser& argParser) {
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

bool D2Cuda::initialize() {
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
    // Everything went on well.
    return true;
}

/** A helper method to check for GPU errors
    Sample usage:
        gpuErrCheck(cudaMalloc(...));
        gpuErrCheck(cudaPeekAtLastError());
*/
#define gpuErrCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
        fprintf(stderr,"GPUassert: %s at %s: %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

/** Method to encode a character to its numerical representation.
    Since this method is only used in buildWordTable, it's annotated with __device__,
    meaning it can only be called on GPU.
*/
__device__ void encode(char c, int &v) {
    if (c == 'A' || c == 'a') v = 0;
    else if (c == 'C' || c == 'c') v = 1;
    else if (c == 'G' || c == 'g') v = 2;
    else if (c == 'T' || c == 't') v = 3;
    else v = -1;
}

/** Method to fill a word table with a sequence's words' numerical representations.
    Since this method is called from host code, it's annotated with __global__.
*/
__global__ void buildWordTable(int *wordTable, char *seq, int numWords, int wordSize, int bitMask) {
    int start = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = start; i < numWords; i += stride) {
        int hash = 0, value;
        bool ignoreHash = false;
        for (int j = i; j < i + wordSize; j++) {
            encode(seq[j], value);
            if (value == -1) {
                ignoreHash = true;
                break;
            }
            hash = ((hash << 2) | value) & bitMask;
        }
        if (!ignoreHash) {
            wordTable[i] = hash;
        }
    }
}

int D2Cuda::setReferenceEST(const EST* est) {
    ASSERT ( est != NULL );
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(est);
    }
    refEST = est;
    // Init ref-est word table
    const std::string s1Str = est->getSequenceString();
    // Calculate length and number of words in sequence
    int s1Length = s1Str.size();
    int numWordsInS1 = s1Length - wordSize + 1;
    numWindowsInS1 = s1Length - frameSize + 1;
    // Allocate the word table on CUDA
    cudaMalloc(&s1WordTable, sizeof(int) * numWordsInS1);
    // Allocate the sequence on CUDA so we can read it inside kernel
    char *s1;
    cudaMalloc(&s1, sizeof(char) * s1Length);
    cudaMemcpy(s1, s1Str.c_str(), sizeof(char) * s1Length, cudaMemcpyHostToDevice);
    // Calculate number of blocks, rounded up
    int numBlocks = (numWordsInS1 + numThreads - 1) / numThreads;
    // Call the actual buildWordTable kernel
    buildWordTable<<<numBlocks, numThreads>>>(s1WordTable, s1, numWordsInS1, wordSize, bitMask);
    // Wait for CUDA to finish and free the sequence memory
    cudaDeviceSynchronize();
    cudaFree(s1);

    return 0;
}

float D2Cuda::getMetric(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    VALIDATE({
        if (otherEST->getID() == refEST->getID()) {
            return 0; // distance to self will be 0
        }
    });
    // OK. Run the actual d2 algorithm
    return (float) runD2(otherEST);
}

/** Method to calculate min D2 score. For each thread running this method, it picks
    a pair of windows from two sequences and calculate D2 score for this pair.
    A delta array is used to keep track of words count differences between two windows.
    Its size is fixed to 4096 since a word largest possible numerical representation is 4095.
    Since this method is called from host code, it's annotated with __global__.
*/
__global__ void getScore(int* s1WordTable, int* s2WordTable, int* minScore, bool* done, int numWindowsInS1, int numWindowsInS2, int numWordsInWindow, int threshold) {
    int start = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int delta[4096];
    // Array values seem to be already initialized with 0
    // Writing to array in GPU is very slow
    // for (int i = 0; i < 4096; i++) {
    //     delta[i] = 0;
    // }
    for (int i = start; i < numWindowsInS1 * numWindowsInS2; i += stride) {
        if (*done) {
            return;
        }
        int s1WindowIndex = i % numWindowsInS1, s2WindowIndex = i / numWindowsInS1, j;
        for (j = s1WindowIndex; j < s1WindowIndex + numWordsInWindow; j++) {
            delta[s1WordTable[j]] += 1;
        }
        for (j = s2WindowIndex; j < s2WindowIndex + numWordsInWindow; j++) {
            delta[s2WordTable[j]] -= 1;
        }
        int score = 0;
        for (j = 0; j < 4096; j++) {
            if (delta[j] != 0) {
                score += delta[j] * delta[j];
                delta[j] = 0;
            }
        }
        atomicMin(minScore, score);
        if (score <= threshold) {
            *done = true;
        }
    }
}

float D2Cuda::runD2(const EST *otherEST) {
    const std::string s2Str = otherEST->getSequenceString();
    // Calculate length and number of words in sequence
    int s2Length = s2Str.size();
    int numWordsInS2 = s2Length - wordSize + 1;
    numWindowsInS2 = s2Length - frameSize + 1;
    // Allocate the word table on CUDA
    cudaMalloc(&s2WordTable, sizeof(int) * numWordsInS2);
    // Allocate the sequence on CUDA so we can read it inside kernel
    char *s2;
    cudaMalloc(&s2, sizeof(char) * s2Length);
    cudaMemcpy(s2, s2Str.c_str(), sizeof(char) * s2Length, cudaMemcpyHostToDevice);
    // Calculate number of blocks, rounded up
    int numBlocks = (numWordsInS2 + numThreads - 1) / numThreads;
    // Call the actual kernel
    buildWordTable<<<numBlocks, numThreads>>>(s2WordTable, s2, numWordsInS2, wordSize, bitMask);
    // Wait for CUDA to finish and free the sequence memory
    cudaDeviceSynchronize();
    cudaFree(s2);

    bool* done;
    int* minScore;
    cudaMalloc(&done, sizeof(bool));
    cudaMallocManaged(&minScore, sizeof(int));
    *minScore = INT_MAX;

    numBlocks = (numWindowsInS1 * numWindowsInS2 + numThreads - 1) / numThreads;
    getScore<<<numBlocks, numThreads>>>(s1WordTable, s2WordTable, minScore, done, numWindowsInS1, numWindowsInS2, numWordsInWindow, threshold);

    cudaDeviceSynchronize();

    int score = *minScore;

    cudaFree(s1WordTable);
    cudaFree(s2WordTable);
    cudaFree(minScore);
    cudaFree(done);

    return (float) score;
}

bool D2Cuda::getAlignmentData(int &alignmentData) {
    // Simply copy the alignment metric that was computed by the last
    // successful call to the analyze() method.
    alignmentData = alignmentMetric;
    // Let the caller know that the alignment data is available.
    return true;
}

#endif
