#ifndef LENGTH_FILTER_CPP
#define LENGTH_FILTER_CPP

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

#include "LengthFilter.h"
#include "ClusterMaker.h"
#include "EST.h"

// Define the static parameters
int LengthFilter::minESTLen = 50;

// The set of arguments for this class.
arg_parser::arg_record LengthFilter::argsList[] = {
    {"--minESTLen", "Minimum EST length to be permitted by this filter",
     &LengthFilter::minESTLen, arg_parser::INTEGER},    
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

LengthFilter::LengthFilter(ClusterMaker *clusterMaker) :
    Filter("lengthFilter", clusterMaker) {
    // Initialize instance variables.
    clusterID = -1;
}


void
LengthFilter::showArguments(std::ostream& os) {
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(LengthFilter::argsList);
    os << ap;
}

bool
LengthFilter::parseArguments(int& argc, char **argv) {
    // Let's process parameters now.
    arg_parser ap(LengthFilter::argsList);
    ap.check_args(argc, argv, false);
    // Ensure values are valid.
    if (minESTLen < 1) {
        std::cerr << filterName
                  << ": Minimum EST length is too small "
                  << "(use --minESTLen option)\n";
        return false;
    }
    // Everything went well
    return true;
}

int
LengthFilter::initialize() {
    ASSERT ( clusterMaker != NULL );
    // Add a dummy cluster for this filter with a suitable name.
    clusterID = clusterMaker->addDummyCluster("Cluster with short ESTs "
                                              "(filtered out by LengthFilter)");
    return 0;
}

int
LengthFilter::runFilter(const int estIndex) {
    ASSERT ( clusterID != -1 );
    ASSERT ((estIndex >= 0) && (estIndex < EST::getESTCount()));
    // Check the length of the specified EST
    const char *seq = EST::getEST(estIndex)->getSequence();
    if ((int) strlen(seq) < minESTLen) {
        // This EST is too short. Filter it out.
        return clusterID;
    }
    // Let this EST through the filter
    return -1;
}

#endif
