#ifndef BATON_ASSEMBLER_MANAGER_CPP
#define BATON_ASSEMBLER_MANAGER_CPP

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

#include "BatonAssemblerManager.h"
#include "ContigMaker.h"
#include "MPIStats.h"
#include "EST.h"

BatonAssemblerManager::BatonAssemblerManager() 
    : BatonAssembler() {
    // Nothing else to be done here for now.
}

BatonAssemblerManager::~BatonAssemblerManager() {
    // Nothing else to be done here for now.
}

bool
BatonAssemblerManager::initialize() {
    // First let the base class do its tasks...
    if (!BatonAssembler::initialize()) {
        // Base class initialization failed
        return false;
    }
    // Everything went fine.
    return true;
}

int
BatonAssemblerManager::assemble() {
    // Cut initial progress entry.
    int estsProcessed = getESTsProcessed();
    reportProgress(estsProcessed);
    // Keep forming contigs until all fragments are processed.
    int refESTidx = -1;
    do {
        // Create a contig using a suitable un-processed fragment as
        // the reference. The following method call to the base class
        // returns a suitable reference EST to be used by the manger
        // and all the worker(s) in one call.
        refESTidx = getReferenceEST();
        if (refESTidx == -1) {
            // All fragments have been processed and the manager can
            // now quit.
            break;
        }
        // Mark reference EST as having been processed.
        EST* const refEST = estList->get(refESTidx, true);
        refEST->setProcessed(true);
        // Create initial contig with reference EST as seed.        
        ContigMaker contig(refESTidx, estList, true);
        // Do some of the alignment on the manager (all alignments
        // happen on the manager if there are zero workers).
        estsProcessed += localAssembly(contig, true);
        // When control drops here, we have a valid consensus sequence
        // to deal with.
        // contig.prettyPrint(std::cout);
        // std::cout << contig.getConsensus() << std::endl;
        // Update progress made thus far
        reportProgress(estsProcessed);
    } while (refESTidx != -1);
    // Processing completed successfully
    return 0;
}

#endif
