#ifndef FILTER_FACTORY_CPP
#define FILTER_FACTORY_CPP

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

#include "FilterFactory.h"
#include "LengthFilter.h"
#include "LCFilter.h"
#include "arg_parser.h"

void
FilterFactory::displayList(std::ostream &os) {
    // Create dummy command-line args to make display prettier and
    // easier.
    arg_parser::arg_record dummy_args[] = {
        {"lengthFilter", "Filter out short ESTs in the FASTA file",
         NULL, arg_parser::STRING},
        {"lcFilter", "Filter out entries with Low Complexity regions",
         NULL, arg_parser::STRING},
        {NULL, NULL, NULL, arg_parser::BOOLEAN}
    };
    arg_parser dummyParser(dummy_args);
    os << dummyParser;
}

Filter*
FilterFactory::create(const char* name, ClusterMaker *clusterMaker) {
    if (name == NULL) {
        return NULL;
    }
    if (clusterMaker == NULL) {
        std::cerr << "A valid cluster maker has  not been specified.\n"
                  << "Use --clusterMaker command line option." << std::endl;
        return NULL;
    }
    
    if (!strcmp("lengthFilter", name)) {
        return new LengthFilter(clusterMaker);
    } else if (!strcmp("lcFilter", name)) {
        return new LCFilter(clusterMaker);
    }
        
    // invalid filter name!
    std::cerr << "Invalid filter name '" << name << "'." << std::endl;
    return NULL;
}

#endif
