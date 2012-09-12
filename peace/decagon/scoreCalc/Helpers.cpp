#ifndef HELPERS_CPP
#define HELPERS_CPP

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

#include "Helpers.h"
#include <cstdio>
#include <cerrno>
#include <sys/stat.h>
#include <iostream>

bool
Helpers::runCmd(const std::string& cmdLine, const int expectedExitStatus) {
    FILE *pipe = NULL;
    // Add redirection to suppress standard error.
    const std::string redirCmdLine = cmdLine + " 2> /dev/null";
    if ((pipe = popen(redirCmdLine.c_str(), "r")) == NULL) {
        // The command was not really found and did not complete
        // successfully.
        return false;
    }
    // Read and ignore outputs from the process (if any)
    char buffer[2048];
    while (fread(buffer, sizeof(char), sizeof(buffer), pipe) != 0);
    // Obtain the exist code from the process.
    const int exitCode = pclose(pipe);
    std::cout << "exitCode[" << cmdLine << "]=" << exitCode << std::endl;
    // Return if the command was successfully run.
    return (exitCode == expectedExitStatus);
}

#endif
