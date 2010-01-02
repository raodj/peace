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
#include "Utilities.h"

// The static instance variables for command line arguments.


// The common set of arguments for all EST analyzers
arg_parser::arg_record ClusterMaker::commonArgsList[] = {
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

ClusterMaker::ClusterMaker(const std::string& myName, ESTAnalyzer *myAnalyzer,
                           const int estIdx, const std::string& outputFile) 
    : name(myName), refESTidx(estIdx), outputFileName(outputFile),
      analyzer(myAnalyzer) {
    // Nothing else to be done for now.
}

ClusterMaker::~ClusterMaker() {}

void
ClusterMaker::showArguments(std::ostream& os) {
    if (commonArgsList[0].arg_text == NULL) {
        // No valid arguments are currently used for all
        // ClusterMakers.  So skip displaying any arguments.
        return;
    }
    
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(ClusterMaker::commonArgsList);
    os << "Common options for all cluster makers are:\n";
    os << ap;
}

bool
ClusterMaker::parseArguments(int& argc, char **argv) {
    arg_parser ap(ClusterMaker::commonArgsList);
    // Process the arguments
    ap.check_args(argc, argv, false);
    // Everything went well.
    return true;
}

// Dummy operator= to keep VC'09 happy.
ClusterMaker& 
ClusterMaker::operator=(const ClusterMaker&) {
    return *this;
}

#endif
