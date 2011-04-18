#ifndef FILTERING_SUB_SYSTEM_CPP
#define FILTERING_SUB_SYSTEM_CPP

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

#include "FilteringSubSystem.h"
#include "FilterFactory.h"

FilteringSubSystem::FilteringSubSystem() {
    filterChain.setSubSystem(this);
}

FilteringSubSystem::~FilteringSubSystem() {
    // Nothing to be done (at least for now)
}

void
FilteringSubSystem::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord Arguments[] = {
        {"", "\nValid arguments for the filtering sub-system are:",
         NULL, ArgParser::INFO_MESSAGE},
        {"--filters", "Names of the filters to use, in order (null for none)",
         &filterNames, ArgParser::STRING_LIST},
        {"", "Valid filter names are:", NULL, ArgParser::INFO_MESSAGE},
        {"", "lengthFilter (Filter out short cDNA entries)",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "lcFilter (Filter out fragments with Low Complexity regions)",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add the valid arguments to the argument parser
    argParser.addValidArguments(Arguments);
    // Have the components add their own parameters.
    filterChain.addCommandLineArguments(argParser);
}

int
FilteringSubSystem::initializeSubSystem(ArgParser& argParser) {
    // Check and setup default entries if none have been specified.
    if (filterNames.empty()) {
        filterNames.push_back("lengthFilter");
        filterNames.push_back("lcFilter");
    }
    // Process one filter name at a time.
    for(size_t filIdx = 0; (filIdx < filterNames.size()); filIdx++) {
        if (filterNames[filIdx] == "null") {
            // A dummy filter name to be ignored.
            break;
        }
        // Create filter and add it to the chain
        Filter *filter =
            FilterFactory::create(filterNames[filIdx], runtimeContext);
        if (filter == NULL) {
            // Break out and return null, which will cause an error
            // and the usage message to show.
            return 1;
        }
        // Have the filter add its own parameters to the list.
        filter->addCommandLineArguments(argParser);
        // Add the filter to our chain.
        filterChain.addFilter(filter);
    }
    // Everything went well.
    return 0;
}

int
FilteringSubSystem::run() {
    // Initialize all the filters.
    if (!filterChain.initialize()) {
        // Initialization failed.
        return 1;
    }    
    // Now let the filter chain run and do filtering.
    if (filterChain.run() != 0) {
        // Error running filter chain.
        return 2;
    }
    // Finally, clear out any dummy entries that may have been added
    // by the filters...
    filterChain.finalize();
    // Everything went on well
    return NO_ERROR;
}

#endif
