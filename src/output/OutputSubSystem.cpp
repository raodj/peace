#ifndef OUTPUT_SUB_SYSTEM_CPP
#define OUTPUT_SUB_SYSTEM_CPP

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

#include "OutputSubSystem.h"
#include "ArgParser.h"
#include "Utilities.h"
#include "RuntimeContext.h"
#include "ClusterMaker.h"
#include "MPIHelper.h"

OutputSubSystem::OutputSubSystem() {
    stdOutput.setSubSystem(this);
    mstWriter.setSubSystem(this);
    clusterWriter.setSubSystem(this);
    samWriter.setSubSystem(this);
    aceWriter.setSubSystem(this);
}

OutputSubSystem::~OutputSubSystem() {
    // Nothing to be done (atleast for now)
}

void
OutputSubSystem::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"", "\nValid arguments for output sub-system are:",
         NULL, ArgParser::INFO_MESSAGE},
        {"--filter-fail", "FASTA file to which ESTs that failed the filter(s) "
         "must be written", &filterFailFile, ArgParser::STRING},
        {"--filter-pass", "FASTA file to which ESTs that pass the filter(s) "
         "are to be written", &filterPassFile, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(Arguments);
    stdOutput.addCommandLineArguments(argParser);
    mstWriter.addCommandLineArguments(argParser);
    clusterWriter.addCommandLineArguments(argParser);
    samWriter.addCommandLineArguments(argParser);
    aceWriter.addCommandLineArguments(argParser);
}

int
OutputSubSystem::initializeSubSystem(ArgParser& UNREFERENCED_PARAMETER(argParser)) {
    if (!stdOutput.initialize()) {
        // Error during standard output stream redirection.
        return 1;
    }
    // Everything went on well
    return NO_ERROR;    
}

int
OutputSubSystem::initializeSubComponents() {
    if (!samWriter.initialize()) {
        // Error during SAM file initialization.
        return 2;
    }
    if (!aceWriter.initialize()) {
        // Error during ACE file initialization.
        return 3;
    }
    // Everything went on well
    return NO_ERROR;
}

void
OutputSubSystem::generateOutputs(const bool UNREFERENCED_PARAMETER(success)) {
    ASSERT( runtimeContext != NULL );

    // The following outputs are generated only by the manager
    if (MPI_GET_RANK() == MANAGER_RANK) {
        const ClusterMaker *clsMaker = runtimeContext->getClusterMaker();
        if ((clsMaker != NULL) && (clsMaker->getMST() != NULL)) {
            mstWriter.write(*clsMaker->getMST());
        }
        
        // Let the cluster writer do its thing.
        clusterWriter.write(runtimeContext);
    }
    // The SAM file writing is performed by all the parallel
    // processes.  It decides if it has necessary information to
    // write the contigs out. So we don't do any extra checks here.
    samWriter.write(runtimeContext);
    // The ACE file writing is performed by all the parallel
    // processes.  It decides if it has necessary information to write
    // the contigs out. So we don't do any extra checks here.
    aceWriter.write(runtimeContext);
    
}

void 
OutputSubSystem::finalizeSubSystem(const bool UNREFERENCED_PARAMETER(success)){
    stdOutput.finalize();
}

#endif
