#ifndef MULTI_WORD_HEURISTIC_CPP
#define MULTI_WORD_HEURISTIC_CPP

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

#include "MultiWordHeuristic.h"
#include "RuntimeContext.h"
#include "HeuristicChain.h"
#include "ArgParser.h"
#include "SubSystem.h"
#include "MPIHelper.h"
#include "ESTCodec.h"
#include "ESTList.h"
#include "EST.h"

// Since the following variables are passed in a template parameter,
// it is defined to be static (to ensure that it has external linkage
// as per the ISO/ANSI standards requirement).
int MultiWordHeuristic::BitMask     = 0;
// Instance variable to track number of bits to shift for building
// hash values.  This is set to 2*(v-1) in initialize method.
int MultiWordHeuristic::bitsToShift = 0;

MultiWordHeuristic::MultiWordHeuristic(HeuristicChain* chain,
				       const std::string& name)
    : Heuristic(name, chain) {
    // Setup default values for parameters
    wordCount = 4;
    wordLen   = 8;
}

MultiWordHeuristic::~MultiWordHeuristic() {
    // Nothing to be done here for now.
}

void
MultiWordHeuristic::addCommandLineArguments(ArgParser& argParser) {
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--mw_wordCount", "Number of matching words in similar fragments",
         &wordCount, ArgParser::INTEGER},
        {"--mw_wordLen", "The number of nucleotides in a word",
         &wordLen, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    // add valid arguments to the command-line parser
    argParser.addValidArguments(ArgsList);
}

bool
MultiWordHeuristic::initialize() {
    // Ensure parameter values are valid.
    if ((wordLen < 1) || (wordLen > 16)) {
        std::cerr << heuristicName
                  << ": Word length must be in the range 1 to 16 (inclusive)"
                  << "(use --mw_wordLen option)\n";
        return false;
    }
    if (wordCount < 1) {
	std::cerr << heuristicName
                  << ": Word count must be greater than 0 (zero)."
                  << "(use --mw_wordCount option)\n";
        return false;
    }
    // Setup the bits to be shifted to create hash values.
    bitsToShift = 2 * (wordLen - 1);
    // Everything went well
    return true;
}

int
MultiWordHeuristic::setReferenceEST(const EST* est) {
    if ((refEST != NULL) && (est != NULL) &&
        (refEST->getID() == est->getID())) {
        // This EST is already the reference EST. Nothing further to
        // do. So just get out.
        return 0;
    }
    // Ensure est is valid.
    if (est == NULL) {
        // Invalid est index.
        return 1;
    }
        
    return 1;
}

void
MultiWordHeuristic::updateWordMap(const EST* est) {
    ASSERT( est != NULL );
    ensurePopulated(est);
    // The bit mask has been defined to be a static to enable passing
    // it as a parameter to the templatized NormalEncoder.
    BitMask = (1 << (wordLen * 2)) - 1;
    // Setup the look up hash table for the reference est.
    refEST = est;
    // Obtain EST sequence for quick access.
    const char* s1 = est->getSequence();    
    // First compute the hash for a single word using a suitable
    // generator.
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    ESTCodec& codec = ESTCodec::getCodec();
    
    int hash        = 0;
    int ignoreMask  = 0;
    const int estID = est->getID();
    for(int i = 0; (i < wordLen - 1); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
    }
    // Fill in the word map. But first setup the table to quickly
    // generate reverse complement values.
    codec.setRevCompTable(wordLen);
    const int End = (int) est->getSequenceLength() - wordLen;
    for (int i = wordLen - 1; (i <= End); i++) {
        hash = encoder(hash, s1[i], ignoreMask);
        if (!ignoreMask) {
            // The ignoreMask is zero implying that this hash does not
            // contain any 'n' bases that are to be ignored. So go
            // ahead and Setup the word map for the regular word
            wordMap[hash][estID] = true;
            // Setup the word map entry for the reverse-complement word.
	    // wordMap[codec.encode2rc(hash)][-estID] = true;
        }
    }
}

bool
MultiWordHeuristic::runHeuristic(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    // Extra sanity checks on uncommon scenarios.
    VALIDATE({
        if (otherEST->getID() == refEST->getID()) {
            return true; // will end up with distance 0, or max similarity
        }
    });
    // When control drops here that means number of matches were not met
    return false;
}

int
MultiWordHeuristic::run() {
    ASSERT ( subSystem != NULL );
    ASSERT ( subSystem->getContext() != NULL );
    ASSERT ( subSystem->getContext()->getESTList() != NULL );    
    // Get convenience reference to global list of ESTs
    ESTList& estList = *(subSystem->getContext()->getESTList());
    // Apply filters on the subset of ESTs that we are responsible for
    // on this process.
    std::cout << "Starting heuristic...\n";
    int startIndex, endIndex;
    getOwnedESTidx(estList.size(), startIndex, endIndex);
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        // Use helper method to update the words associated with the
        // filters.
        EST *est = estList[estIdx];
        if (!est->hasBeenProcessed()) {
	    // Process the read.
	    updateWordMap(est);
        }
    }
    std::cout << "End heuristic...\n";
    // Now participate in global iterative broadcast chain. Possibly
    // this loop can be replaced by a MPI all-to-all broadcast that
    // expense of increased memory footprint.
    if (MPI_GET_SIZE() > 1) {
        // ClusterMaker *clusterMaker =
	// subSystem->getContext()->getClusterMaker();
        // allToAllBroadcast(clusterMaker);
    }
    // All operations proceeded successfully
    return 1;
}

#endif
