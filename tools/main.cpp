#ifndef MAIN_CPP
#define MAIN_CPP

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "ArgParser.h"
#include "ShowAlignment.h"
#include "PrimeFeatures.h"
#include "DrawMatrix.h"
#include "ShowMST.h"
#include "Matcher.h"
#include "Align.h"

// The list of valid tools and associated usage.
const char* ToolInfo[] = {
    "ShowAlignment","Show graphical view of EST alignments",
    "DrawMatrix",   "Draw graphical view of D2  matrix of scores",
    "ShowMST",      "Translate MST to a bracketed form for rendering",
    "Align",        "Simple alignment generation tool using MST alignment data",
    "Matcher",      "Score contigs against source DNA/Transcripts",
    "primes",       "Print n-dimensional features based on prime numbers",
    // Add new tool names and description before this line.
    NULL, NULL
};

/** \func showUsage

    A simple helper method to show usage message and the list of
    tools associated with searums tools.

    \note This method is meant to be in local scope rather than global
    scope.
*/
static void showUsage() {
    std::cout << "PEACE Tools Software Suite.\n"
              << "Usage: pTool --tool <toolName>\n"
              << "where valid tool names are:\n";
    for(int idx = 0; (ToolInfo[idx] != NULL); idx += 2) {
        std::cout << "\t" << ToolInfo[idx] << ": "
                  << ToolInfo[idx + 1] << std::endl;
    }
}

/** \func main

    <p> The global main method for the PEACE tools software suite. The
    main method performs the task of identifying the specific tool to
    run and pass control to the specific class.  Accordingly, it only
    consumes the "--tool" option which is typically the first option
    specified.  However, it is not required that --tool be the first
    option. </p>

    <p>Once a specific PEACE tool hass be identified instantiated the
    main method simply passes control over to that specific tool's
    main() method.</p>
*/
int
main(int argc, char* argv[]) {
    std::string tool;
    // Create the list of valid arguments to be used by the arg_parser.
    ArgParser::ArgRecord arg_list[] = {
        {"--tool", "Name of the tool to be run",
         &tool, ArgParser::STRING},
        {"", "", NULL, ArgParser::INVALID}
    };

    // Get the argument parser to parse and consume the global
    // options.  Based on the options supplied, various variables will
    // be set to appropriate values.
    ArgParser ap(arg_list);
    ap.parseArguments(argc, argv, false);
    // Ensure a valid tool name has been specified.
    int toolCode = 0;
    for(toolCode = 0; (ToolInfo[toolCode] != NULL); toolCode += 2) {
        if (ToolInfo[toolCode] == tool) {
            // Found a match!
            break;
        }
    }
    
    if (ToolInfo[toolCode] == NULL) {
        // A valid tool was not specified.
        showUsage();
        // Return nothing.
        return 1;
    }

    // Based on the tool chain to other tools.
    int retVal = 0;
    switch(toolCode) {
    case 0: retVal = ShowAlignment::main(argc, argv);
        break;
    case 2: retVal = DrawMatrix::main(argc, argv);
        break;
    case 4: retVal = ShowMST::main(argc, argv);
        break;
    case 6: retVal = Align::main(argc, argv);
        break;
    case 8: retVal = Matcher::main(argc, argv);
        break;
    case 10: retVal = PrimeFeatures::main(argc, argv);
        break;
    default:
        std::cout << "Invalid or unhandled tool name specified.\n";
    }
    
    // return the result back to the caller.
    return retVal;
}

#endif
