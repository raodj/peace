#ifndef BATON_ASSEMBLER_CPP
#define BATON_ASSEMBLER_CPP

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

#include "BatonAssembler.h"
#include "ArgParser.h"
#include "MPIHelper.h"
#include "MPIStats.h"
#include "ArgParser.h"
#include "ContigMaker.h"
#include "RuntimeContext.h"
#include "BatonAlignmentInfo.h"

#include <sstream>

BatonAssembler::BatonAssembler() : Assembler("baton") {
    // Nothing else to be done here.
}

BatonAssembler::~BatonAssembler() {
    // Nothing to be done here.
}

void
BatonAssembler::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add its arguments first
    Assembler::addCommandLineArguments(argParser);
    // Let the cache adds its parameters.
    blCache.addCommandLineArguments(argParser);
    // Let the sequence aligner add its parameters
    aligner.addCommandLineArguments(argParser);
}

bool
BatonAssembler::initialize() {
    // First let the base class load the cDNA data
    if (!Assembler::initialize()) {
        // Base class initialization failed. This is a serious
        // problem.
        return false;
    }
    // If the runtime context has an baton analyzer, we will report
    // a warning as we don't really use the analyzer.
    ASSERT( getContext() != NULL );
    if (getContext()->getAnalyzer() != NULL) {
        std::cerr << "Warning: BatonAssembler does not use an analyzer.\n"
                  << "Use the command line option '--analyzer null' to "
                  << "avoid this warning.\n";
    }
    aligner.setBLCache(&blCache);
    // Setup cross references in our personal aligner
    aligner.setSubSystem(subSystem);
    // Now initialize the aligner so it can do its validation(s).
    if (!aligner.initialize()) {
        // Analyzer initialization failed. This is a serious problem.
        return false;
    }
    // Initialization finished successfully.
    return true;
}

void
BatonAssembler::finalize() {
    // Display MPI usage statistics.
    std::cout << "MPI Statistics from rank " << MPI_GET_RANK() << std::endl;
    MPIStats::displayStats(std::cout);    
}

int
BatonAssembler::localAssembly(ContigMaker& startContig, const bool isManager) {
    // Track number of ESTs processed to report progress information
    // in manager process.
    int estsProcessed = 0;
    if ((estsProcessed = growContig(startContig)) == 0) {
        // No additional ESTs were found to match. Nothing further to
        // be done.
        return 0;
    }
    // Shortcut to the window size for convenience.
    const int WindowSize = aligner.getWindowSize();
    // Next, repeatedly grow initialContig to the left until no new
    // ESTs are added to the sub contigs. Each sub-contig created is
    // merged into the startContig
    while (startContig.hasLeftOverhang()) {
        // Create a new contig maker for exploring to the left
        ContigMaker leftContig =
            startContig.getOverhangContigMaker(true, WindowSize, isManager);
        // Now check unprocessed, locally-owned ESTs for similarity to
        // this contig and align highly-similar ESTs
        estsProcessed += growContig(leftContig);
        // Merge new contig with the main contig (causing main contig
        // to grow to left)
        startContig.merge(leftContig);
    }

    // Next check and grow the initialContig to the right and merge
    // the new contigs into startContig
    while (startContig.hasRightOverhang()) {
        // Create a new contig maker for exploring to the left
        ContigMaker rightContig =
            startContig.getOverhangContigMaker(false, WindowSize, isManager);
        // Grow the right contig
        estsProcessed += growContig(rightContig);
        // Merge with the main contig
        startContig.merge(rightContig);
    }
    // Now that the contig creation is complete, we create and add the
    // contig.
    createAddContig(startContig, isManager);
    // Return number of ESTs processed
    return estsProcessed;
}

int
BatonAssembler::growContig(ContigMaker& contigMaker) {
    // Get the subset of fragments to be analyzed and assembled by
    // this process. First setup the reference fragment for assembly
    if (contigMaker.getReferenceEST() != NULL) {
        aligner.setReferenceEST(contigMaker.getReferenceEST());
    } else {
        aligner.setReferenceEST(contigMaker.getReferenceSeq());
    }
    // Get the locally-owned ESTs that are under the purview of this
    // process.
    int startIndex, endIndex;
    getOwnedESTidx(estList->size(), startIndex, endIndex);
    // Check and process each unprocessed fragment, while tracking
    // number of ESTs added to the contig maker.
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        EST* const otherEST = estList->get(estIdx, true);
        if (otherEST->hasBeenProcessed()) {
            // This EST entry can be ignored as it has already been
            // processed.
            continue;
        }
        // Check and align estIdx and reference sequence, if they are
        // similar. The alignment information is obtained in the
        // BatonAlignmentInfo object defined below.
        BatonAlignmentInfo alignInfo;
        if (aligner.align(otherEST, alignInfo)) {
            // An alignment was found!  Add this fragment to the
            // contig to let contig maker grow.
            contigMaker.add(alignInfo);
            // Flag EST as having been processed.
            otherEST->setProcessed(true);
        }
    }
    // Form global contig by merging information from parallel
    // processes
    const int estsAdded = contigMaker.formGlobalContig();
    // Return number of ESTs added.
    return estsAdded;
}

int
BatonAssembler::getReferenceEST() {
    // Get the subset of fragments under the purview of this process
    // to locate the longest, local, un-processed fragment.
    int startIndex, endIndex;
    getOwnedESTidx(estList->size(), startIndex, endIndex);
    // Check and process each unprocessed fragment.
    int maxLenIdx = -1;  // Index of longest EST
    int maxLen = -1;  // The longest EST thusfar
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        EST* const otherEST = estList->get(estIdx);
        if (otherEST->hasBeenProcessed()) {
            // This EST entry can be ignored as it has already been
            // processed.
            continue;
        }
        // Track the longest est and its index
        if (maxLen < otherEST->getSequenceLength()) {
            // Found a longer fragment to use
            maxLen    = otherEST->getSequenceLength();
            maxLenIdx = estIdx;
        }
    }
    // Now distribute the information to all processes while selecting
    // the largest value.
    int localMaxData[2]  = {maxLen, maxLenIdx};
    int globalMaxData[2] = {maxLen, maxLenIdx};  // Convenience in no-MPI-case
    MPI_ALL_REDUCE(localMaxData, globalMaxData, 1, MPI_TYPE_2INT, MPI_MAXLOC);
    // Return the global, longest, available, EST index as the initial
    // reference cDNA fragment for contig formation.
    return globalMaxData[1];
}

void
BatonAssembler::createAddContig(const ContigMaker& contigMaker,
                                const bool isManager) {
    // A static counter to generate unique conting numbers.
    static int contigCounter = 1;
    // Add an empty contig to the list. The newly added conting will
    // be populated further below.
    contigList.push_back(Contig());
    Contig& contig = contigList.back();
    // Create a unique identifier for this contig.
    std::ostringstream contigID;
    contigID << "Contig" << contigCounter;
    contig.setId(contigID.str());
    // Increase contig counter
    contigCounter++;    
    // Next, get the contig maker to pouplate the contig information
    // us.
    contigMaker.populateContig(contig, isManager, isManager);
    // Finally, notify all interested listeners that a new contig has
    // been made.
    const bool removeContig = notifyContigListeners(contig, false);
    if (removeContig) {
        // The contig has been successfully processed. No need to hold
        // onto it or its data. First, we unpopulate cDNA information
        // to minimize memory footprint.
        const AlignmentInfoList& aiList = contig.getAlignments();
        for(size_t i = 0; (i < aiList.size()); i++) {
            EST *est = estList->get(aiList[i].getESTIndex());
            ASSERT ( est != NULL );
            est->unpopulate();
            // Unpouplate baton lists for this EST from our cache as
            // we no longer need it.
            blCache.unpopulate(est);
        }
        contigList.pop_back();
    }
}

int
BatonAssembler::getESTsProcessed(const int procID) const {
    // First gather the number of ESTs already processed (due filters etc.)
    int localESTsProcessed = estList->getProcessedESTCount();
    // The following contain global counts at procID after MPI_Reduce
    // call below (at other processes the value remains unchanged).
    int estsProcessed      = localESTsProcessed;
    MPI_REDUCE(&localESTsProcessed, &estsProcessed, 1, MPI_TYPE_INT,
               MPI_OP_SUM, procID);
    // Return the counts.
    return estsProcessed;
}

#endif
