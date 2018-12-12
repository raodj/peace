#ifndef MAIN_CPP
#define MAIN_CPP

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

#include "InteractiveConsole.h"
#include "framework/ArgParser.h"
#include "framework/MPIHelper.h"
#include "PEACE.h"


/** \fn int main(int argc, char *argv[])

    \brief The main function that launches all the activites.

    <p> The global main method for the ESTAnalyzer software. The main
    method performs the task of identifying and running PEACE either
    in interactive or non-interactive mode.  Accordingly, it only
    consumes the "--interactive" option which is typically the first
    option specified.  However, it is not required that --analyzer be
    the first option. </p>

    \param[in] argc The number command-line arguments specified to the
    executable.  This variable indicates the number of entries in argv
    parameter.

    \param[in] argv The array of command-line arguments.
*/
int
main(int argc, char* argv[]) {
    // Perform any mpi initialization as needed
    MPI_INIT(&argc, &argv);
    
    // Values of the following variables are processed by the argument
    // parser further below.
    bool interactive = false;
    // Create the list of valid arguments to be used by the arg_parser.
    const ArgParser::ArgRecord argList[] = {
        {"--interactive", "Run PEACE in interactive mode.",
         &interactive, ArgParser::BOOLEAN},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Get the argument parser to parse and consume just the
    // --interactive flag (for now).
    ArgParser ap(argList);
    ap.parseArguments(argc, argv, false);

    // Now create the PEACE object for this run
    PEACE peace;
    int exitCode = 0;
    if (interactive) {
        // Launch PEACE interactive console.
        InteractiveConsole console(&peace);
        console.processCommands(argc, argv);
    } else {
        // Run PEACE in non-interactive mode.
        exitCode = peace.run(argc, argv, false);
    }
    // Windup the PEACE object.
    peace.finalize(false);
    // Shutdown MPI.
    MPI_FINALIZE();
    // Return result of analysis (usually 0 to indicate no errors).
    return exitCode;
}

#endif
