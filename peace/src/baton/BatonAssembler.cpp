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

#include "ESTAnalyzerFactory.h"
#include "BatonAssembler.h"
#include "BatonAnalyzer.h"
#include "MPIHelper.h"
#include "MPIStats.h"

// Static member variables used to define the arg records below.
char* BatonAssembler::progFileName;

// The set of arguments specific for baton assembler.  Several
// parameters are handled by the BatonAnalyzer class.
arg_parser::arg_record BatonAssembler::argsList[] = {
    {"--progress", "Log assembly progress in a file (used by GUI)",
     &BatonAssembler::progFileName, arg_parser::STRING},            
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};


BatonAssembler::BatonAssembler(const std::string& outputFileName)
    : Assembler("baton", outputFileName),
      analyzer(dynamic_cast<BatonAnalyzer*>(ESTAnalyzerFactory::create("baton", 0, outputFileName))) {
    // Ensure our baton-analyzer has been successfully created.
    ASSERT( analyzer != NULL );
}

BatonAssembler::~BatonAssembler() {
    if (analyzer != NULL) {
        // Delete the analyzer as we no longer need it.
        delete analyzer;
    }
}

void
BatonAssembler::showArguments(std::ostream& os) {
    // Let base class show its arguments first
    Assembler::showArguments(os);
    // Let the baton analyzer that we use by default show its
    // arguments next
    ASSERT ( analyzer != NULL );
    analyzer->showArguments(os);
    // Finally, show arguments that are specific to this assembler.
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(BatonAssembler::argsList);
    os << ap;
}

bool
BatonAssembler::parseArguments(int& argc, char **argv) {
    // First, let base class process common arguments first
    if (!Assembler::parseArguments(argc, argv)) {
        // Error processing common arguments. Bail out
        return false;
    }
    // Let the baton analyzer that we use by default process all of
    // its arguments next
    ASSERT ( analyzer != NULL );
    if (!analyzer->parseArguments(argc, argv, true)) {
        // Error processing arguments for the baton analyzer.
    }
    // Finally, process the arguments that are specific to this
    // assembler via a suitable arg_parser.
    arg_parser ap(BatonAssembler::argsList);
    ap.check_args(argc, argv, false);
    // Check if values of command line parameters are within
    // acceptable bounds.

    // Everything went well
    return true;
}

int
BatonAssembler::initialize() {
    // First let the base class load the cDNA data
    int retVal = Assembler::initialize();
    // If the data was successfully loaded, initialize the analyzer
    if (retVal == 0) {
        ASSERT ( analyzer != NULL );
        retVal = analyzer->initialize(false);
    }
    if (retVal != 0) {
        // Return overall error back to the caller
        return retVal;
    }
    // Everything went well.
    return 0;
}

void
BatonAssembler::getOwnedESTidx(int& startIndex, int& endIndex) {
    // Now compute the subset of ESTs to be processed by this
    // process. The following code (until the for-loop) evenly
    // distributes the ESTs between all the processes.
    const int ProcessCount   = MPI_GET_SIZE();
    const int ESTsPerProcess = EST::getESTList().size() / ProcessCount;
    const int ExtraESTs      = EST::getESTList().size() % ProcessCount;
    const int MyRank         = MPI_GET_RANK();
    // First figure out the starting and ending EST this processs is
    // responsible for further use.
    startIndex = MyRank * ESTsPerProcess;
    // Evenly divide up any extra ESTs between the workers.
    if (MyRank <= ExtraESTs) {
        // The previous processes have one extra ESTs as number of
        // ESTs are not evenly divisible by number of processes.  So
        // account for this.
        startIndex = ((ESTsPerProcess + 1) * MyRank);
    } else {
        startIndex += ExtraESTs;
    }
    // Compute the last est index this process owns.
    endIndex = startIndex + ESTsPerProcess;
    if (MyRank < ExtraESTs) {
        // This process owns one extra EST as ESTs are not evenly
        // divisible by the number of processes.
        endIndex++;
    }
}

void
BatonAssembler::localAssembly(const int refESTidx,
                              AlignmentInfoList& alignList) {
    // Get the subset of fragments to be analyzed and assembled by
    // this process. First setup the reference fragment for assembly
    ASSERT( analyzer != NULL );
    analyzer->setReferenceEST(refESTidx);
    int startIndex, endIndex;
    getOwnedESTidx(startIndex, endIndex);
    // Check and process each unprocessed fragment.
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        if (EST::getEST(estIdx)->hasBeenProcessed()) {
            // This EST entry can be ignored as it has already been
            // processed.
            continue;
        }
        // Check and align estIdx and refESTidx, if they are
        // similar. The alignment information is obtained in the
        // AlignmentInfo object defined below.
        AlignmentInfo alignInfo;
        if (analyzer->align(estIdx, alignInfo)) {
            // An alignment was found, then the information is added
            // Save the alignment for future reference
            alignList.push_back(alignInfo);
            // Flag EST as having been processed.
            EST::getEST(estIdx)->setProcessed(true);
        }
    }
}

#endif
