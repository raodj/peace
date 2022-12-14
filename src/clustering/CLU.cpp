#ifndef CLU_CPP
#define CLU_CPP

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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include "CLU.h"
#include "ResultLog.h"
#include "ESTList.h"
#include "ArgParser.h"

#include <cmath>
#include <queue>
#include <algorithm>
#include <cstring>

CLU::CLU() : FWAnalyzer("clu") {
    abundanceFraction = 10; // 10 is default from CLU implementation.    
    // Initialize the array CharToInt for mapping A, T, C, and G to 0,
    // 1, 2, and 3 respectively.
    memset(CharToInt, 0, sizeof(CharToInt));
    // Initialize values for specific base pairs.
    CharToInt[(int) 'A'] = CharToInt[(int) 'a'] = 0;
    CharToInt[(int) 'G'] = CharToInt[(int) 'g'] = 1;
    CharToInt[(int) 'C'] = CharToInt[(int) 'c'] = 2;
    CharToInt[(int) 'T'] = CharToInt[(int) 't'] = 3;
}

CLU::~CLU() {
    // Clearing out EST list also frees up the custom CLUESTData
    // associated with the EST's thereby freeing up any reference and
    // complement hash maps.
}

void
CLU::addCommandLineArguments(ArgParser& argParser) {
    // Let base class initialize itself.
    FWAnalyzer::addCommandLineArguments(argParser);
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--abdFrac", "Abundance Fraction (default=10)",
         &abundanceFraction, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Use a arg parser object to conveniently display common options.
    argParser.addValidArguments(ArgsList);
}

bool
CLU::initialize() {
    if (!FWAnalyzer::initialize()) {
        // Error occured when initializing base class.  This is no good.
        return false;
    }
    // Ensure abundance fraction is valid.
    if ((abundanceFraction < 1) || (abundanceFraction > 100)) {
        std::cerr << getName()
                  << ": Abundance fraction must be in the range 0 to 100"
                  << "(use --abdFrac option)\n";
        return false;
    }
    
    // Initialize the custom CLU data to hold the hash maps for each
    // EST.  Note that the hash maps by themselves are not yet
    // populated.
    for(int id = 0; (id < estList->size()); id++) {
        estList->get(id, false)->setCustomData(new CLUESTData());
    }
    // Everything went on well.
    return true;
}

void
CLU::dumpEST(ResultLog& log, const EST* est, const bool isReference) {
    if (isReference) {
        const char* start=(htmlLog ? "<font face=\"courier\" color=red>" : "");
        const char* end  =(htmlLog ? "</font>" : "");
        log.report(est->getInfo(), "%s%s%s", "n/a",
                   start, est->getSequence(), end);
    } else {
        const char* start = (htmlLog ? "<font face=\"courier\">" : "");
        const char* end   = (htmlLog ? "</font>" : "");
        log.report(est->getInfo(), "%s%s%s", "%f",
                   start, est->getSequence(), end, est->getSimilarity());
    }
}

int
CLU::setReferenceEST(const EST *est) {
    if (est == NULL)  {
        // Invalid est provided.  This is an error condition.
        return 1;
    }
    // Save the reference ESTindex for future use.    
    refEST = est;
    // Everything went well.
    return 0;
}

void
CLU::buildHashMaps(const EST *est) {
    ASSERT ( est != NULL );
    ASSERT ( est->getCustomData().get() != NULL );
    // Check if the tables have already been created.
    CLUESTData *customData =
        dynamic_cast<CLUESTData*>(est->getCustomData().get());
    ASSERT ( customData != NULL );
    if (customData->referenceHashMap != NULL) {
        // Hash maps have already been created for this EST.  Nothing
        // further to be done.
        return;
    }
    // Obtain the base pair sequence for this EST.
    std::string sequence = est->getSequence();
    ASSERT ( sequence.size() > 0 );
    
    // Create the reference hash map.
    createCLUHashMap(customData->referenceHashMap, sequence.c_str());
    // Now reverse the sequence in prep to create complement sequence.
    std::reverse(sequence.begin(), sequence.end());
    const char Complements[] = {'T', 'C', 'G', 'A'};
    const int SeqLen = (int) sequence.size();
    for(int i = 0; (i < SeqLen); i++) {
        sequence[i] = Complements[(int) CharToInt[(int) sequence[i]]];
    }
    // Create the complement hash table.
    createCLUHashMap(customData->complementHashMap, sequence.c_str());
}

void
CLU::createCLUHashMap(int* &hashMap, const char *sequence) {
    // Create a reference hash map if needed.
    const int MapSize = (int) pow(4.0, wordSize);
    if (hashMap == NULL) {
        hashMap = new int[MapSize];
    }
    // Clear out any old entries in the reference hash map.
    memset(hashMap, 0, sizeof(int) * MapSize);
    // First compute the hash for a single word.
    ASSERT ( sequence != NULL );
    int hash = 0;
    int i;
    for(i = 0; (i < wordSize); i++) {
        hash <<= 2;
        hash  += CharToInt[(int) sequence[i]];
    }
    // Setup the hash entry for the first word.
    hashMap[hash]++;
    // Now that the hash for the first word has been computed, the
    // hash for subsequent words can be easily derived from it using
    // bit-wise operations as the word simply slides acorss the
    // sequence.
    const int SeqLen  = (int) strlen(sequence) - wordSize;
    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    const int BitMask = (1 << (wordSize * 2)) - 1;
    for (; (i < SeqLen); i++) {
        hash <<= 2;
        hash  += CharToInt[(int) sequence[i]];
        hash  &= BitMask;
        hashMap[hash]++;
    }
    
    // Clean up abundant, non-informative, and low-complexity oligo
    // sequences from the reference hash map using a helper method.
    filterHashMap(hashMap, SeqLen);
}

void
CLU::filterHashMap(int *hashMap, const int sequenceLength) {
    // First, zero out all simple oligos from consideration.  The
    // simple oligos are sequences of the form:
    // "AAAAAACCCCCCGGGGGGTTTTTT".  Each of these low-complexity oligos
    // are zeroed out in the reference hash map.  To ease hash
    // generation the low-complexity oligo sequence is first generated.
    
    // Note that the first set of AA..A is ignored as its hash is
    // always 0 (zero) as value of each 'A' is 0!
    char *LowCompSeq = reinterpret_cast<char*>(alloca(wordSize * 3));
    memset(LowCompSeq + wordSize * 0, 1, wordSize);
    memset(LowCompSeq + wordSize * 1, 2, wordSize);
    memset(LowCompSeq + wordSize * 2, 3, wordSize);

    // Start with the initial hash value of 0 as the first set of
    // AAA..A is always 0 (zero) as the value of each 'A' is 0!
    int hash         = 0;
    hashMap[hash]    = 0;
    const int SeqLen = wordSize * 3;
    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.    
    const int BitMask= (1 << (wordSize * 2)) - 1;
    for(int i = 0; (i < SeqLen); i++) {
        hash <<= 2;
        hash  += LowCompSeq[i];
        hash  &= BitMask;
        // Zero out the hash hashMap entry for this low complexity oligo
        hashMap[hash] = 0;
    }
    // Zero-out abundant sequences as they are not very informative.
    const int MapSize         = (int) pow(4.0, wordSize);
    const int AbundanceFactor = sequenceLength / abundanceFraction;
    for(int i = 0; (i < MapSize); i++) {
        if (!hashMap[i]) {
            // If the entry is zero there is nothing further to do.
            continue;
        }
        if (hashMap[i] > AbundanceFactor) {
            hashMap[i] = 0;
        } else {
            hashMap[i] = 1;
        }
    }
}

int
CLU::getSimilarity(const int* const hashMap,
                   const char* const sequence) const {
    // Create a categorized distribution for computing similarity.
    int *cd = reinterpret_cast<int*>(alloca(sizeof(int) * frameSize));
    // Initialize entries in categorized distribution (cd) to 0.
    memset(cd, 0, sizeof(int) * frameSize);
    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    const int BitMask = (1 << (wordSize * 2)) - 1;
    // Set up a pointer to iterate over base pair sequences.
    const char* bp = sequence;

    // Compute raw hash scores for first word...
    int hash = 0;
    for(int i = 0; (i < wordSize); i++) {
        hash <<= 2;
        hash  += CharToInt[(int) *bp++];
    }
    // Save the hash in a queue to rapidly adjust hash values as the
    // frame slides along the given sequence.
    std::queue<int> prevHashValues;
    prevHashValues.push(hash);
    // Compute the running categorized distribution entry index for
    // the first frame of the given sequence.
    unsigned int cdIndex = hashMap[hash];
    for(int i = 1; (i < frameSize); i++) {
        hash   <<= 2;
        hash    += CharToInt[(int) *bp++];
        hash    &= BitMask;
        cdIndex += hashMap[hash];
        // Save hash value to adjust running values later on...
        prevHashValues.push(hash);
    }
    ASSERT ( (int) prevHashValues.size() == frameSize );
    // Now using the hash and cdIndex values for first frame compute
    // the values for other frames by suitably adjusting the values.
    const int SeqLen = (int) strlen(sequence) - frameSize;
    for(int i = frameSize; (i < SeqLen); i++) {
        // Ensure cdIndex never exceeds the number of cateogries.
        cdIndex = (cdIndex >= (unsigned) frameSize) ? (frameSize - 1) : cdIndex;
        // Track number of frames in a given category..
        cd[cdIndex]++;
        // Using hash for previous frame compute hash values for next
        // frame.
        hash <<= 2;
        hash  += CharToInt[(int) *bp++];
        hash  &= BitMask;
        // Save hash value to adjust running cdIndex values later on...
        prevHashValues.push(hash);
        // Adjucst running cdIndex value next.
        cdIndex += hashMap[hash];
        cdIndex -= hashMap[prevHashValues.front()];
        // Remove the entry already used from the list.
        prevHashValues.pop();
        // Enusre prevListSize remains sane.
        ASSERT ( (int) prevHashValues.size() == frameSize );
    }

    // The CLU implementation does this last step.  It is unclear why
    // this extra step is necessary.  But this code is faithfully
    // replicating CLU behavior (bug-for-bug if this is the case).    
    // Ensure cdIndex never exceeds the number of cateogries.
    cdIndex = (cdIndex >= (unsigned) frameSize) ? (frameSize - 1) : cdIndex;
    // Track number of frames in a given category..
    cd[cdIndex]++;
    
    // Finally, multiply each category distribution value with a
    // pre-calculated weight and sum the products to obtain the final
    // similarity metric.
    const int Weight[] = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                           0,  0,  1,  3,  7, 15, 27, 34, 52, 47};
    ASSERT ( sizeof(Weight) / sizeof(int) == frameSize);
    int similarity = 0;
    for(int i = 0; (i < frameSize); i++) {
        similarity += (cd[i] * Weight[i]);
    }
    // Returns similarity back to the caller.
    return similarity;
}

float
CLU::getMetric(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    if (otherEST->getID() == refEST->getID()) {
        // Comparing reference to itself makes no sense in CLU.
        // Simply return 0 (zero) in this case.
        return 0;
    }
    // Obtain the two EST's to be compared.
    const EST *estSq = refEST;
    const EST *estSb = otherEST;
    // As per CLU's strategy always use the longer sequence as Sq
    if (strlen(estSb->getSequence()) < strlen(estSq->getSequence())) {
        // Sb has longer sequence than Sq.  So swap flip them.
        const EST *temp = estSq;
        estSq     = estSb;
        estSb     = temp;
    }
    // First build the reference table for the query EST, estSq
    buildHashMaps(estSq);
    // Check to ensure the tables have been created and we have ready
    // access to it.
    CLUESTData *customData =
        dynamic_cast<CLUESTData*>(estSq->getCustomData().get());
    ASSERT ( customData != NULL );
    ASSERT ( customData->referenceHashMap  != NULL );
    ASSERT ( customData->complementHashMap != NULL );
    // Now obtain reference and complement similarity metrics
    const char* otherSeq = estSb->getSequence();
    ASSERT ( otherSeq != NULL );
    const int refSim  = getSimilarity(customData->referenceHashMap, otherSeq);
    const int compSim = getSimilarity(customData->complementHashMap, otherSeq);
    // Return the best of refSim or compSim.
    return (float) ((refSim > compSim) ? refSim : compSim);
}

CLU::CLUESTData::~CLUESTData() {
    if (referenceHashMap != NULL) {
        delete [] referenceHashMap;
    }
    if (complementHashMap != NULL) {
        delete [] complementHashMap;
    }
}

#endif
