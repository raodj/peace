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

#include "ScoreCalculator.h"
#include "ACEReader.h"
#include <iostream>

/** \fn int main(int argc, char *argv[])

    \brief The main function that launches all the activites.

    <p> The global main method for the Score Calculator. The main
    function merely delegates control to an instnace of the
    ScoreCalculator class by calling ScoreCalculator::run() with
    necessary parameters. </p>

    \param[in] argc The number command-line arguments specified to the
    executable.  This variable indicates the number of entries in argv
    parameter.

    \param[in] argv The array of command-line arguments.
*/
int
main(int argc, char* argv[]) {
    // Delegate all the tasks to a score calculator object.
    ScoreCalculator sc;
    try {
        // Let the score calculator do its job
        sc.run(argc, argv);
    } catch (const DecagonException& de) {
        // Catch application-specific exceptions.
        std::cerr << de;
        return de.getErrorCode();
    }
    // Everything went on fine
    return 0;
}

#endif
