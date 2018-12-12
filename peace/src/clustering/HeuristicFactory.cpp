#ifndef HEURISTIC_FACTORY_CPP
#define HEURISTIC_FACTORY_CPP

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

#include "HeuristicFactory.h"
#include "UVSampleHeuristic.h"
#include "TVHeuristic.h"
#include "OldTVHeuristic.h"
#include "MultiWordHeuristic.h"
#include "PrimesHeuristic.h"
#include "ArgParser.h"

void
HeuristicFactory::addCommandLineInfo(ArgParser& argParser) {
    // Create dummy command-line args to make display prettier and
    // easier.
    const ArgParser::ArgRecord DummyArgs[] = {
        {"", "mw : The Multi Word heuristic",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "uv : The u/v sample heuristic",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "tv : The t/v heuristic (u/v is automatically included)",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "oldTV: A deprecated version of t/v heuristic",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "primes: A prime number based heuristic",
         NULL, ArgParser::INFO_MESSAGE},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(DummyArgs);
}

Heuristic*
HeuristicFactory::create(const std::string& name,  HeuristicChain *chain) {
    if (name.empty()) {
        return NULL;
    }
    if (name == "uv") {
        return new UVSampleHeuristic(chain);
    } else if (name == "tv") {
        return new TVHeuristic(chain);
    } else if (name == "oldTV") {
        return new OldTVHeuristic(chain);
    } else if (name == "mw") {
	return new MultiWordHeuristic(chain);
    } else if (name == "primes") {
        return new PrimesHeuristic(chain);
    }
    
    // invalid heuristic name!
    std::cerr << "Invalid heuristic name '" << name << "'." << std::endl;
    return NULL;
}

#endif
