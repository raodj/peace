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
#include "RuntimeContext.h"
#include "ArgParser.h"

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
        return false;
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
BatonAssembler::localAssembly(const EST* refEST,
                              AlignmentInfoList& alignList) {
    // Get the subset of fragments to be analyzed and assembled by
    // this process. First setup the reference fragment for assembly
    aligner.setReferenceEST(refEST);
    int startIndex, endIndex;
    getOwnedESTidx(estList->size(), startIndex, endIndex);
    // Check and process each unprocessed fragment.
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        EST* const otherEST = estList->get(estIdx, true);
        if (otherEST->hasBeenProcessed()) {
            // This EST entry can be ignored as it has already been
            // processed.
            continue;
        }
        // Check and align estIdx and refESTidx, if they are
        // similar. The alignment information is obtained in the
        // AlignmentInfo object defined below.
        AlignmentInfo alignInfo;
        if (aligner.align(otherEST, alignInfo)) {
            // An alignment was found, then the information is added
            // Save the alignment for future reference
            alignList.push_back(alignInfo);
            // Flag EST as having been processed.
            otherEST->setProcessed(true);
        }
    }
}

#endif
