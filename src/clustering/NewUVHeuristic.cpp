#ifndef NEW_UV_HEURISTIC_CPP
#define NEW_UV_HEURISTIC_CPP

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

#include "NewUVHeuristic.h"
#include "HeuristicChain.h"
#include "ESTCodec.h"
#include "ArgParser.h"
#include "EST.h"

// Since the following two variables are passed in a template parameter,
// it is defined to be static (to ensure that it has external linkage
// as per the ISO/ANSI standards requirement).
int NewUVHeuristic::BitMask           = 0;
// Instance variable to track number of bits to shift for building
// hash values.  This is set to 2*(v-1) in initialize method.
int NewUVHeuristic::bitsToShift = 0;

NewUVHeuristic::NewUVHeuristic(HeuristicChain* chain, const std::string& name)
    : Heuristic(name, chain) {
    // Initialize hash table arrays
    s1WordMap   = NULL;
    s1RCWordMap = NULL;

    // Set defaults for the dynamic parameters
    v           = 8;
    u           = 6;
    wordShift   = 8;
    passes      = 3;
    refESTLen   = 0;
}

NewUVHeuristic::~NewUVHeuristic() {
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
NewUVHeuristic::addCommandLineArguments(ArgParser& argParser) {
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--uv_v", "v (length of common words)",
         &v, ArgParser::INTEGER},
        {"--uv_u", "u (number of v-word matches for similar fragments)",
         &u, ArgParser::INTEGER},
        {"--uv_wordShift", "Words to shift/skip to make heuristic faster",
         &wordShift, ArgParser::INTEGER},
        {"--uv_passes", "Number of Passes (default=3)",
         &passes, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };    
    // Use a arg parser object to conveniently display common options.
    argParser.addValidArguments(ArgsList);
}

bool
NewUVHeuristic::initialize() {
    // Validate the parameters.
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
    } else if (passes < 1) {
        std::cerr << heuristicName
                  << ": passes must be >= 1"
                  << "(use --uv_passes option)\n";
        return false;
    }
    // Setup the standard data structures
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
NewUVHeuristic::setReferenceEST(const EST* estS1) {
    ASSERT ( estS1 != NULL );
    refEST = estS1;
    // Ensure initialize() method has done its job
    ASSERT( s1WordMap   != NULL );
    ASSERT( s1RCWordMap != NULL );
    
    // The bit mask has been defined to be a static to enable passing
    // it as a parameter to the templatized NormalEncoder.
    BitMask = (1 << (v * 2)) - 1;
    // Setup the look up hash table for the reference est.
    // Initialize s1 word map for the new reference EST
    const int MapSize = (1 << (v * 2));

    // Initialize word maps to false throughout.
    memset(s1WordMap,   0, sizeof(char) * MapSize);
    memset(s1RCWordMap, 0, sizeof(char) * MapSize);
    // Obtain EST sequence for quick access.
    ensurePopulated(estS1);
    const char* s1 = estS1->getSequence();
    refESTLen = estS1->getSequenceLength();
    
    // First compute the hash for a single word using a suitable
    // generator.
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    ESTCodec& codec = ESTCodec::getCodec();
    
    int hash  = 0, ignoreMask = 0;
    for(int i = 0; (i < v - 1); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
    }
    // Fill in the word map. But first setup the table to quickly
    // generate reverse complement values.
    codec.setRevCompTable(v);
    const int End = estS1->getSequenceLength();
    for (int i = v - 1; (i <= End); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
        if (!ignoreMask) {
            // If ignoreMask is zero, it tells us that the hash has no
            // values associated with base pairs tagged as 'n'. So use
            // it to setup the word map for the regular word
            s1WordMap[hash] = 1;
            // Setup the word map entry for the reverse-complement word.
            s1RCWordMap[codec.encode2rc(hash)] = 1;
        }
    }
    return 0; // everything went well
}

void
NewUVHeuristic::computeHash(const EST* estS2) {
    // Pre-compute hash values that will be used for comparisons
    // to reduce duplicate computations.
    ASSERT ( estS2  != NULL );
    ensurePopulated(estS2);
    const char *sq2  = estS2->getSequence();
    ASSERT ( (int) strlen(sq2) == estS2->getSequenceLength() );
    ASSERT ( sq2 != NULL );
    const int End    = estS2->getSequenceLength() - v;
    
    // Get the codec for encoding/decoding operations
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    // Obtain the vector from the hash_map and ensure it has
    // sufficient space to hold all the hash values.
    std::vector<unsigned short>& hashList = uvCache[estS2->getID()];
    hashList.reserve(End / wordShift + 1);
    
    // Now, go through the EST sq2 and build hash values
    for (register int start = 0; (start < End); start += wordShift) {
        // Compute the hash for the next v words.
        unsigned short hash = 0;
        int ignoreMask      = 0;
        const int endBP     = start + v;
        for(register int i  = start; (i < endBP); i++) {
            hash = (unsigned short) encoder(hash, sq2[i], ignoreMask);
        }
        // Store the hash value if usable
        if (!ignoreMask) {
            // If ignoreMask is zero it tells us that the hash does
            // not have 'n' in it and is good to be used.
            hashList.push_back(hash);
        }
    }
}

bool
NewUVHeuristic::runHeuristic(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    // Extra sanity checks on uncommon scenarios.
    VALIDATE({
        if (otherEST->getID() == refEST->getID()) {
            // First clear the hint for the cluster maker
            HeuristicChain::getHeuristicChain()->setHint(HeuristicChain::MST_RC,
                                                         0);
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
    const int hashSize = otherHash.size();
    if (hashSize == 0) {
        // No valid words, therefore this pair need not be analyzed further.
        return false;
    }
    // go through the otherHash and track number of matching words
    // Initialize local variables.
    register int numMatches = 0, numRCmatches = 0;
    register int factor = u / passes;
    for (register int pass = 0; (pass < passes); pass++) {
        for (register int word = pass; (word < hashSize); word+=passes) {
            // Track number of positive and reverse-complement matches on
            // this hash word
            numMatches   += s1WordMap  [otherHash[word]];
            numRCmatches += s1RCWordMap[otherHash[word]];
        }
        if (pass != passes-1) {
            register int numMatchesReq = factor*(pass+1);
            if ((numMatches < numMatchesReq) && (numRCmatches < numMatchesReq)) {
                // Not enough matches, so break out immediately
                // printHists(refEST, otherEST, bestMatchIsRC, false);
                return false;
            }
        }
    }

    // Print some information for analysis purposes.
    // printf("uv(%d, %d) = %d, %d\n", refEST->getID(), 
    //       otherEST->getID(), numMatches, numRCmatches);
    
    // Check if we had sufficient matches to warrant further
    // operations
    if ((numMatches >= u) || (numRCmatches >= u)) {
        // Set the flag to indicate if the normal or the reverse
        // complement version of checks yielded the best result only
        // if needed.
        bestMatchIsRC = (numMatches < numRCmatches);
        // Setup a hint for D2.
        heuristicChain->setHint(HeuristicChain::D2_DO_RC, bestMatchIsRC);
        // Setup a hint for the MST Cluster Maker (-1 = RC, 1 = no RC)
        heuristicChain->setHint(HeuristicChain::MST_RC, bestMatchIsRC ? -1 : 1);
        // printHists(refEST, otherEST, bestMatchIsRC, true);        
        // return success indication
        return true;
    }
    // When control drops here that means number of matches were not met
    // printHists(refEST, otherEST, bestMatchIsRC, false);    
    return false;
}

    
#include <sstream>
#include <fstream>

int shift  = 0, mask = 0;

void
NewUVHeuristic::genHist(const EST* est1, const int nMers,
                        std::vector<int>& hist,
                        std::vector<int>& rcHist) {
    const int buckets = (1 << (nMers * 2));
    hist.resize(buckets);
    rcHist.resize(buckets);
    for(int i = 0; (i < buckets); i++) {
        hist[i]   = 0;
        rcHist[i] = 0;
    }
    // First compute the hash for a single word using a suitable
    // generator.
    shift = 2 * (nMers - 1);
    mask  = (1 << (nMers * 2)) - 1;
    ESTCodec::NormalEncoder<shift, mask> encoder;
    ESTCodec& codec = ESTCodec::getCodec();
    // Initialize the first window
    int hash  = 0, ignoreMask = 0;
    const char* s1 = est1->getSequence();
    for(int i = 0; (i < nMers - 1); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
    }
    // Fill in the word map. But first setup the table to quickly
    // generate reverse complement values.
    codec.setRevCompTable(nMers);
    const int End = est1->getSequenceLength();
    for (int i = nMers - 1; (i <= End); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
        // Track the number of entries
        hist[hash]++;
        // Setup the word map entry for the reverse-complement word.
        rcHist[codec.encode2rc(hash)]++;
    }
}

void
NewUVHeuristic::printHist(const EST* est1,
                          const EST* est2,
                          const bool revCompMatch,
                          const bool result, int nMers) {
    std::ostringstream fileName;
    fileName << "est" << est1->getID() << "_est" << est2->getID() << "_"
             << (revCompMatch ? "RC" : "NN")
             << (result ? "_Pass_" : "_Fail_")
             << "Nmer" << nMers << ".dat";
    std::ofstream dataFile(fileName.str().c_str());
    
    std::vector<int> est1Hist, est1RcHist, est2Hist, est2RcHist;
    genHist(est1, nMers, est1Hist, est1RcHist);
    genHist(est2, nMers, est2Hist, est2RcHist);
    for(size_t i = 0; (i < est1Hist.size()); i++) {
        dataFile << i
                 << "\t" << est1Hist[i] << "\t" << est1RcHist[i]
                 << "\t" << est2Hist[i] << "\t" << est2RcHist[i] << std::endl;
    }
    dataFile.close();
    std::cout << "Histogram written to " << fileName.str() << std::endl;
    std::cout << "Press ENTER to continue..." << std::flush;
    std::string dummy;
    std::getline(std::cin, dummy);
    std::getline(std::cin, dummy);
}

void
NewUVHeuristic::checkHist(const EST* est1,
                          const EST* est2,
                          const bool revCompMatch,
                          const bool result, int nMers) {
    
    std::vector<int> est1Hist, est1RcHist, est2Hist, est2RcHist;
    genHist(est1, nMers, est1Hist, est1RcHist);
    genHist(est2, nMers, est2Hist, est2RcHist);
    // Search in histogram for matching entries with just 3 entries.
    int matches = 0;
    for(size_t i = 0; (i < est1Hist.size()); i++) {
        if (((est1Hist[i] > 4) || (est1RcHist[i] > 4)) &&
            (est2Hist[i] > 4)) {
            // This pair is matched!
            matches++;
        }
    }
    if (result != (matches > 0)) {
        std::cout << "** Disagreement EST#" << est1->getID()
                  << " <-> EST#" << est2->getID()
                  << " matches = " << matches
                  << " and result = "
                  << (result ? "true" : "false")
                  << ", revCompMatch = "
                  << (revCompMatch ? "true" : "false") << std::endl;
        if (result) {
            std::cout << "View histogram (y/n)? " << std::flush;
            std::string ans;
            std::cin  >> ans;
            if (ans == "y") {
                printHist(est1, est2, revCompMatch, result, nMers);
            }
        }
    }
}

void
NewUVHeuristic::printHists(const EST* est1,
                           const EST* est2,
                           const bool revCompMatch,
                           const bool result) {
    checkHist(est1, est2, revCompMatch, result, 6);
}

#endif
