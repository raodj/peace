#ifndef CLUSTER_MAKER_FACTORY_CPP
#define CLUSTER_MAKER_FACTORY_CPP

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

#include "ClusterMakerFactory.h"
#include "ArgParser.h"
#include "AdaptiveMSTClusterMaker.h"
#include "NonAdaptiveMSTClusterMaker.h"
#include "TransMSTClusterMaker.h"

void
ClusterMakerFactory::addCommandLineInfo(ArgParser& argParser) {
    // Create dummy command-line args to make display prettier and
    // easier.
    const ArgParser::ArgRecord DummyArgs[] = {
        {"", "mst   : Adaptive MST-based Cluster Maker (best option)",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "amst  : Adaptive MST-based Cluster Maker (aka mst)",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "na-mst: Non-adaptive MST-based Cluster Maker for Sanger data",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "tmst  : MST-based Cluster Maker with Transitivity",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(DummyArgs);
}

ClusterMaker*
ClusterMakerFactory::create(const std::string& name, ESTAnalyzer *analyzer) {
    if (analyzer == NULL) {
        std::cerr << "A EST analyzer has not been specified.\n"
                  << "Use --analyzer command line option." << std::endl;
        return NULL;
    }
    
    if ((name == "mst") || (name == "amst")) {
        return new AdaptiveMSTClusterMaker(analyzer);
    } else if (name == "tmst") {
        return new TransMSTClusterMaker(analyzer);
    } else if (name == "na-mst") {
        return new NonAdaptiveMSTClusterMaker(analyzer);
    }
    
    // invalid analyzer name!
    std::cerr << "Invalid analyzer name." << std::endl;
    return NULL;
}

#endif
