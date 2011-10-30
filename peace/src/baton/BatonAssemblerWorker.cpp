#ifndef BATON_ASSEMBLER_WORKER_CPP
#define BATON_ASSEMBLER_WORKER_CPP

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

#include "BatonAssemblerWorker.h"
#include "ContigMaker.h"
#include "MPIHelper.h"
#include "MPIStats.h"

BatonAssemblerWorker::BatonAssemblerWorker()
    : BatonAssembler() {
    // Nothing else to be done here for now.
}

BatonAssemblerWorker::~BatonAssemblerWorker() {
    // Nothing else to be done here for now.
}

int
BatonAssemblerWorker::assemble() {
    // The operations in the following loop must be coordinated with
    // the sequence of operations in the manager process as well.
    // The reference EST index that is updated in the loops below.
    ASSERT( estList != NULL );
    // Coordinate with manager to report initial set of ESTs already
    // processed at this worker.
    getESTsProcessed();
    
    int refESTidx = -1;
    do {
        // First get the index of the reference sequence with which
        // the subset of ESTs assigned to this worker must be analyzed
        // and assembled. The following method call to the base class
        // returns a suitable reference EST to be used by the manger
        // and all the worker(s).
        refESTidx = getReferenceEST();
        if (refESTidx == -1) {
            // All fragments have been processed and the worker can
            // now quit.
            break;
        }
        // Mark reference EST as having been processed.
        EST* const refEST = estList->get(refESTidx, true);
        refEST->setProcessed(true);
        // Create initial contig with reference EST as seed.
        ContigMaker contig(refESTidx, estList, false);                
        // Do some of the alignment on the worker using the common
        // base class helper.
        localAssembly(contig, false);
    } while (refESTidx != -1);
    // All done. Return zero indicating success
    return 0;
}

#endif
