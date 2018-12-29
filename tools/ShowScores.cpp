#ifndef SHOW_SCORES_CPP
#define SHOW_SCORES_CPP

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include "ShowScores.h"
#include "Common.h"
#include "EST.h"
#include "ESTList.h"
#include "ESTAnalyzer.h"
#include "ArgParser.h"

int
ShowScores::main(int argc, char *argv[]) {
    int  refRead    = 0;     // Reference EST reads
    bool showOptions= false; // Flag to indicat if help is to be displayed

    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--ref-read", "Reference read for comparison/metric generation",
         &refRead, ArgParser::INTEGER},
        {"--options", "Lists options for this tool",
         &showOptions, ArgParser::BOOLEAN},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    ArgParser ap;
    // Let base class add some defaults for us.
    Tool::addCmdLineArgs("ShowAlignment", ap);
    // Add our custom arguments
    ap.addValidArguments(arg_list);
    // Validate the command-line options.
    ap.parseArguments(argc, argv, false);    
    if (showOptions) {
        std::cout << ap << std::endl;
        // Nothing further to be done.
        return 0;
    }
    // Create the show scores object for us to use
    ShowScores ss(argc, argv);
    ss.showScores(refRead);
    // Everything went well.
    return 0;
}

void
ShowScores::showScores(const int refID) {
    // Get the list of ESTs to be processed
    ESTList& estList = getESTList();
    // Get the analyzer to use from PEACE
    ESTAnalyzer& analyzer = getAnalyzer();
    // Setup the initial reference read
    analyzer.setReferenceEST(estList.get(refID, true));
    // Print the scores for each read
    for (int i = 0; (i < estList.size());i ++) {
        std::cout << refID << ", " << i << ", "
                  << analyzer.analyze(estList.get(i, true))
                  << std::endl;
    }
}

// Construtor
ShowScores::ShowScores(int argc, char* argv[]) : Tool(argc, argv) {
    // Set primes to be used in the base helper class
}

#endif
