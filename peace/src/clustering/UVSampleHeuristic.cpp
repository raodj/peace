#ifndef UV_SAMPLE_HEURISTIC_CPP
#define UV_SAMPLE_HEURISTIC_CPP

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

#include "UVSampleHeuristic.h"
#include "HeuristicChain.h"
#include "ArgParser.h"
#include "ESTCodec.h"
#include "EST.h"

// Since the following variables are passed in a template parameter,
// it is defined to be static (to ensure that it has external linkage
// as per the ISO/ANSI standards requirement).
int UVSampleHeuristic::BitMask     = 0;
// Instance variable to track number of bits to shift for building
// hash values.  This is set to 2*(v-1) in initialize method.
int UVSampleHeuristic::bitsToShift = 0;

UVSampleHeuristic::UVSampleHeuristic(HeuristicChain* chain,
                                     const std::string& name)
    : Heuristic(name, chain) {
    // Initialize hash table arrays
    s1WordMap   = NULL;
    s1RCWordMap = NULL;
    // Setup default values for parameters
    u           = 4;
    v           = 8;
    wordShift   = 16;

}

UVSampleHeuristic::~UVSampleHeuristic() {
    if (s1WordMap != NULL) {
        delete [] s1WordMap;
    }
    if (s1RCWordMap != NULL) {
        delete [] s1RCWordMap;
    }
    // Clear out all the entires in the uvCache as they are no longer
    // needed.
    uvCache.clear();
}

void
UVSampleHeuristic::addCommandLineArguments(ArgParser& argParser) {
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--uv_u", "u (number of v-word matches for similar fragments)",
         &u, ArgParser::INTEGER},
        {"--uv_v", "v (length of common words)",
         &v, ArgParser::INTEGER},
        {"--uv_wordShift", "Words to Shift/skip to make heuristic faster",
         &wordShift, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    // add valid arguments to the command-line parser
    argParser.addValidArguments(ArgsList);
}

bool
UVSampleHeuristic::initialize() {
   // Ensure parameter values are valid.
    if (wordShift < 1) {
        std::cerr << heuristicName
                  << ": Word shift must be >= 1"
                  << "(use --uv_wordShift option)\n";
        return false;
    }
    if ((u < 1) || (v < 1)) {
        std::cerr << heuristicName
                  << ": u and v must be >= 1"
                  << "(use --uv_u and --uv_v options)\n";
        return false;
    }
    // Perform initialization and setup tasks
    const int MapSize = (1 << (v * 2));
    // Create the normal encoded word map
    s1WordMap = new char[MapSize];
    // Create the reverse-complement encoded word map
    s1RCWordMap = new char[MapSize];
    // Setup the bits to be shifted to create hash values.
    bitsToShift = 2 * (v - 1);
    // Everything went well
    return true;
}

int
UVSampleHeuristic::setReferenceEST(const EST* est) {
    if ((refEST != NULL) && (est != NULL) &&
        (refEST->getID() == est->getID())) {
        // This EST is already the reference EST. Nothing further to
        // do. So just get out.
        return 0;
    }
    
    // The bit mask has been defined to be a static to enable passing
    // it as a parameter to the templatized NormalEncoder.
    BitMask = (1 << (v * 2)) - 1;
    // Ensure est is valid.
    if (est == NULL) {
        // Invalid est index.
        return 1;
    }
    // Setup the look up hash table for the reference est.
    refEST = est;
    // Initialize s1 word map for the new reference EST
    const int MapSize = (1 << (v * 2));

    // Initialize word maps to false throughout.
    memset(s1WordMap,   0, sizeof(char) * MapSize);
    memset(s1RCWordMap, 0, sizeof(char) * MapSize);
    // Obtain EST sequence for quick access.
    ensurePopulated(est);
    const char* s1 = est->getSequence();
    
    // First compute the hash for a single word using a suitable
    // generator.
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    ESTCodec& codec = ESTCodec::getCodec();
    
    register int hash = 0;
    int ignoreMask    = 0;
    for(int i = 0; (i < v - 1); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
    }
    // Fill in the word map. But first setup the table to quickly
    // generate reverse complement values.
    codec.setRevCompTable(v);
    const int End = (int) est->getSequenceLength() - v;
    for (int i = v - 1; (i <= End); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
        if (!ignoreMask) {
            // The ignoreMask is zero implying that this hash does not
            // contain any 'n' bases that are to be ignored. So go
            // ahead and Setup the word map for the regular word
            s1WordMap[hash] = 1;
            // Setup the word map entry for the reverse-complement word.
            s1RCWordMap[codec.encode2rc(hash)] = 1;
        }
    }
    return 0; // everything went well
}

void
UVSampleHeuristic::computeHash(const EST* est) {
    // Pre-compute values that will be used in the for-loop below to
    // reduce duplicate computations.
    ASSERT ( est  != NULL );
    ensurePopulated(est);
    const char *sq2  = est->getSequence();
    ASSERT ( sq2 != NULL );
    ASSERT ( (int) strlen(sq2) == est->getSequenceLength() );
    const int End    = std::max(est->getSequenceLength(),
                                (int) est->getSequenceLength() - v);
    ASSERT ( End > 0 );
    
    // Get the codec for encoding/decoding operations
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    // Obtain the vector from the hash_map and ensure it has
    // sufficient space to hold all the hash values.
    std::vector<unsigned short>& hashList = uvCache[est->getID()];
    hashList.reserve(End / wordShift + 1);
    
    // Now, go through the EST sq2 and build hash values
    for (register int start = 0; (start < End); start += wordShift) {
        // Compute the hash for the next v words.
        unsigned short hash = 0;
        int ignoreMask      = 0;
        const int endBP    = start + v;
        for(register int i = start; (i < endBP); i++) {
            hash = (unsigned short) encoder(hash, sq2[i], ignoreMask);
        }
        // Store the hash value if it did not have a 'n' in it.
        if (!ignoreMask) {
            // This hash does not have a 'n' in it. So use it.
            hashList.push_back(hash);
        }
    }
}

bool
UVSampleHeuristic::runHeuristic(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    // Extra sanity checks on uncommon scenarios.
    VALIDATE({
        if (otherEST->getID() == refEST->getID()) {
            return true; // will end up with distance 0, or max similarity
        }
    });
    
    // Get otherEST's hash values from the uvCache. If an entry for
    // otherEST is not present in uvCache then build one.
    UVHashTable::iterator cacheEntry = uvCache.find(otherEST->getID());
    if (cacheEntry == uvCache.end()) {
        // An entry does not exist. Create one and cache it for future
        // references.
        computeHash(otherEST);
        // Update our iterator for use below
        cacheEntry = uvCache.find(otherEST->getID());
    }
    // Obtain a reference to the hash list from the cache.
    const std::vector<unsigned short>& otherHash = cacheEntry->second;
    const int hashSize = (int) otherHash.size();
    if (hashSize == 0) {
        // No valid words, therefore this pair need not be analyzed further.
        return false;
    }
    // go through the otherHash and track number of matching words
    // Initialize local variables.
    register int numMatches = 0, numRCmatches = 0;
    for (register int word = 0; (word < hashSize); word++) {
        // Track number of positive and reverse-complement matches on
        // this hash word
        numMatches   += s1WordMap  [otherHash[word]];
        if (numMatches > 10) {
            numMatches++;
        }
        numRCmatches += s1RCWordMap[otherHash[word]];
    }
    
    // Check if we had sufficient matches to warrant further
    // operations
    if ((numMatches >= u) || (numRCmatches >= u)) {
        // Set the flag to indicate if the normal or the reverse
        // complement version of checks yielded the best result
        bestMatchIsRC = (numMatches < numRCmatches);
        // Setup a hint for D2.
        heuristicChain->setHint(HeuristicChain::D2_DO_RC, bestMatchIsRC);
        // return success indication
        return true;
    }
    // When control drops here that means number of matches were not met
    return false;
}

#endif
