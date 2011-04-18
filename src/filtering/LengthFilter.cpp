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
#include "ArgParser.h"
#include "Utilities.h"
#include "EST.h"

LengthFilter::LengthFilter(RuntimeContext *runtimeContext) :
    Filter("lengthFilter", runtimeContext) {
    // Initialize instance variables.
    clusterID = -1;
    minESTLen = 50;
}

void
LengthFilter::addCommandLineArguments(ArgParser& argParser) {
    // The set of arguments for this filter
    const ArgParser::ArgRecord LengthFilterArgs[] = {
        {"--minESTLen", "Minimum EST length to be permitted by length filter",
         &minESTLen, ArgParser::INTEGER},    
        {"", "", NULL, ArgParser::INVALID}
    };
    // Use a arg parser object to add command line arguments specific
    // for this filter's operations    
    argParser.addValidArguments(LengthFilterArgs);
}

int
LengthFilter::initialize() {
    // Add a dummy cluster for this filter with a suitable name.
    clusterID = addDummyCluster("Cluster with short ESTs "
                                "(filtered out by LengthFilter)");
    return 0;
}

bool
LengthFilter::runFilter(const EST *est, int& clusterID) {
    ASSERT ( est != NULL );
    // Check the length of the specified EST
    const int seqLen = est->getSequenceLength();
    if (seqLen < minESTLen) {
        // This EST is too short. Filter it out.
        clusterID = this->clusterID;
        return true;
    }
    // Let this EST through the filter
    return false;
}

#endif
