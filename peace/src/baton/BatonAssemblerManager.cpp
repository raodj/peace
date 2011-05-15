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
#include "MPIHelper.h"
#include "MPIStats.h"
#include "EST.h"

// Reference to the constant zero to make code more readable
#define MANAGER_RANK 0

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
    // Populate the list of all fragments to be assembled in the
    // estsToProcess vector (used in while-loop below)
    setupESTsToBeProcessed();
    // Cut initial progress entry if requested
    reportProgress(0);
    // The reference EST index that is updated in the loops below.
    int refESTidx = -1;
    // Keep forming contigs until all fragments are processed.
    while (!estsToProcess.empty()) {
        // Create a contig using the first un-processed fragment as
        // the reference.  First create a contig maker using the first
        // available fragment.
        ContigMaker contig(*estsToProcess.begin(), estList);
        estsToProcess.erase(estsToProcess.begin());
        while ((refESTidx = contig.nextReference()) != -1) {
            // Distribute reference index to all workers (zero or more)
            MPI_BCAST(&refESTidx, 1, MPI_INT, MANAGER_RANK);
            // Mark reference EST as having been processed.
            EST* const refEST = estList->get(refESTidx, true);
            refEST->setProcessed(true);
            // Do some of the alignment on the manager (all alignments
            // happen on the manager if there are zero workers).
            AlignmentInfoList alignmentList;
            localAssembly(refEST, alignmentList);
            // Populate the consensus sequence with the information
            // from local assembly.
            if (alignmentList.size() > 0) {
                addToContig(contig, &alignmentList[0], alignmentList.size());
            }
            // Gather information from various workers into the contig
            // object to build consensus sequence.
            getWorkerAlignments(contig);
        }
        // When control drops here, we have a valid consensus sequence
        // to deal with.
        contig.formContig(true);
        // Update progress made thus far
        reportProgress(estList->size() - estsToProcess.size());
    }
    // Finally broadcast -1 as reference index to indicate that all
    // sequences have been assembled.    
    ASSERT ( refESTidx == -1 );
    MPI_BCAST(&refESTidx, 1, MPI_INT, MANAGER_RANK);
    return 0;
}

void
BatonAssemblerManager::addToContig(ContigMaker& contig,
                                   const AlignmentInfo* alignList,
                                   const int listSize) {
    ASSERT ( alignList != NULL );
    for(int index = 0; (index < listSize); index++) {
        // Add alignment information to the overall consensus list.
        contig.add(alignList[index]);
        // Remove the added list from the local set of fragments to be
        // processed.
        estsToProcess.erase(estsToProcess.find(alignList[index].getESTIndex()));
    }
}

void
BatonAssemblerManager::getWorkerAlignments(ContigMaker& contig) {
    const int WorkerCount = MPI_GET_SIZE() - 1;
    for(int workerID = 0; (workerID < WorkerCount); workerID++) {
        MPI_CODE({
                // Get alignment list from any worker that has the
                // alignment list ready for the manager to
                // consume. First determine the size of the alignment
                // list.
                MPI_STATUS msgInfo;
                MPI_PROBE(MPI_ANY_SOURCE, ALIGN_INFO_LIST, msgInfo);
                // Allocate necessary memory to copy the data.  This
                // is a bit of a hack that we are exchanging
                // structures as flat character arrays and then
                // type-casting them to appropriate data type.
                const size_t dataSize = msgInfo.Get_count(MPI_TYPE_CHAR);
                ASSERT ( dataSize >= sizeof(AlignmentInfo) );
                char *buffer = new char[dataSize];
                ASSERT ( buffer != NULL );
                // Receive the alignment info list as an array.
                MPI_RECV(buffer, dataSize, MPI_TYPE_CHAR, msgInfo.Get_source(),
                         ALIGN_INFO_LIST);
                // Type-cast array to suitable object. This works
                // because we assume SPMD and heterogeneous clusters
                const AlignmentInfo* alignList =
                    reinterpret_cast<AlignmentInfo*>(buffer);
                const int listSize = dataSize / sizeof(AlignmentInfo);
                ASSERT ( listSize > 0 );
                ASSERT ( alignList != NULL );
                // let helper method process the list for us while
                // ignoring last dummy entry.
                addToContig(contig, alignList, listSize - 1);
                // Now free the memory for message as we no longer need it
                delete [] buffer;
            });
    }
}

void
BatonAssemblerManager::setupESTsToBeProcessed() {
    // Create the set of ESTs to be processed across the board
    const int ESTCount = estList->size();
    for(int estIdx = 0; (estIdx < ESTCount); estIdx++) {
        if (!estList->get(estIdx)->hasBeenProcessed()) {
            // Insert the index of EST in the list of ESTs to be
            // processed by this process.
            estsToProcess.insert(estIdx);
        }
    }
}

#endif
