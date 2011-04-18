#ifndef RUNTIME_CONTEXT_CPP
#define RUNTIME_CONTEXT_CPP

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

#include "RuntimeContext.h"
#include "Version.h"
#include "Utilities.h"
#include <cstdlib>

RuntimeContext::RuntimeContext() {
    // Add some of the standard/default information to the runtime
    // configuration information.
    std::string peaceInfo(VERSION);
    peaceInfo = peaceInfo + " (Released on: " + RELEASE_DATE + ")";
    addConfig("PEACE Version", peaceInfo);
    // Add time when configuration was created
    char now[128];
    getTime(now);
    addConfig("Runtime Timestamp", now);
    // Reset pointers to shared data structures
    estList      = NULL;
    analyzer     = NULL;
    clusterMaker = NULL;
    assembler    = NULL;
    estListener  = NULL;
}

void
RuntimeContext::addCmdLineInfo(const int argc, char*argv[]) {
    // Create a string for the global arguments to add to
    // configuration information.
    std::string args = "";
    for(int i = 0; (i < argc); i++) {
        args += argv[i];
        args += " ";
    }
    // Add global arguments to the configuration
    addConfig("Command Line", args);
}

void
RuntimeContext::printConfig(std::ostream& os) const {
    StringStringMap::const_iterator curr;
    for(curr = config.begin(); (curr != config.end()); curr++) {
        os << "# " << curr->first << ": " << curr->second << std::endl;
    }
}

void
RuntimeContext::setESTList(ESTList* list) {
    if (estList != NULL) {
        std::cerr << "A valid pointer to the shared list of cDNA fragments "
                  << "has already been set in the RuntimeContext and it "
                  << "cannot be changed. Aborting!" << std::endl;
        abort();
    }
    // Set the reference only once
    ASSERT( list != NULL );
    estList = list;
};

void
RuntimeContext::setAnalyzer(ESTAnalyzer* estAnalyzer) {
    if (analyzer != NULL) {
        std::cerr << "A valid pointer to the shared ESTAnalyzer object "
                  << "has already been set in the RuntimeContext and it "
                  << "cannot be changed. Aborting!" << std::endl;
        abort();
    }
    // Set the reference only once
    ASSERT( estAnalyzer != NULL );
    analyzer = estAnalyzer;
};

void
RuntimeContext::setClusterMaker(ClusterMaker* clusterer) {
    if (clusterMaker != NULL) {
        std::cerr << "A valid pointer to the shared ClusterMaker object "
                  << "has already been set in the RuntimeContext and it "
                  << "cannot be changed. Aborting!" << std::endl;
        abort();
    }
    // Set the reference only once
    ASSERT( clusterer != NULL );
    clusterMaker = clusterer;
};

void
RuntimeContext::setAssembler(Assembler* assembler) {
    if (this->assembler != NULL) {
        std::cerr << "A valid pointer to the shared Assembler object "
                  << "has already been set in the RuntimeContext and it "
                  << "cannot be changed. Aborting!" << std::endl;
        abort();
    }
    // Set the reference only once
    ASSERT( assembler != NULL );
    this->assembler = assembler;
};

void
RuntimeContext::setESTListListener(ESTListListener* listener) {
    estListener = listener;
}

#endif
