#ifndef CLUSTER_MAKER_CPP
#define CLUSTER_MAKER_CPP

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

#include "ClusterMaker.h"
#include "SubSystem.h"
#include "Utilities.h"
#include "RuntimeContext.h"
#include "ESTAnalyzer.h"

ClusterMaker::ClusterMaker(const std::string& myName, ESTAnalyzer *myAnalyzer) 
    : Component(myName), analyzer(myAnalyzer) {
    // Nothing else to be done for now.
}

ClusterMaker::~ClusterMaker() {
    // Nothing else to be done for now.
}

void
ClusterMaker::addCommandLineArguments(ArgParser& argParser) {
    // Currently nothing important to be added to the argument parser
    // but for a message.
    const ArgParser::ArgRecord CommonArgsList[] = {
        {"", "Arguments for " + getName() + " cluster maker:",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(CommonArgsList);
}

bool
ClusterMaker::initialize() {
    if (!Component::initialize()) {
        // Base class inititalization failed
        return false;
    }
    // Setup our convenience pointer to shared ESTList
    ASSERT( subSystem != NULL );
    ASSERT( subSystem->getContext() != NULL );
    ASSERT( subSystem->getContext()->getESTList() != NULL );
    estList = subSystem->getContext()->getESTList();
    // Inititalize the ESTAnalyzer (which inturn inititalizes the
    // heuristic chain and the ReferenceSetManager).
    if (!analyzer->initialize()) {
        // Error occured during initialization. Bail out.
        return false;
    }
    // Return success.
    return true;
}

// Dummy operator= to keep VC'09 happy.
ClusterMaker& 
ClusterMaker::operator=(const ClusterMaker&) {
    return *this;
}

#endif
