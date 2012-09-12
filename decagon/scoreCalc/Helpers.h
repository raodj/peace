#ifndef HELPERS_H
#define HELPERS_H

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

/** \file Helpers.h

	\brief A static class with various utility methods.

	This is a static class that encapsulates various helper/utility
	methods that are used by various classes in Decagon.
*/

#include <string>

/** Class to contain various helper/utility methods.

    This class is a non-instantiable wrapper to contain various
    helper/utility methods.  All the public methods in this class are
    static -- consequently, they can be directly invoked without
    instantiating an object of this class.
*/
class Helpers {
public:
    /** Helper method to run a given command.

        This method essentially attempts to run the specified command
        and verifies that the command runs successfully. This method
        reads all standard outputs from the program and simply
        discards it.

        \param[in] cmdLine The command-line (including the executable
        name) to be run.  If this method is being used to detect if a
        specific executable is available, then the command-line should
        \i typically attempt to perform a quick operation (such as
        displaying help) that can be quickly finished to verify
        existence of the desired executable. However, quick running
        commands are not a requirement and even long running commands,
        if just a success/failure result is sufficient.

        \param[in] expectedExitStatus The expected exit code from the
        process. The default success status code is typically zero.
        
        \return This method returns \c true if the specified command
        executable was found on the path, the command could be
        successfully run, and the exit code matches the specified exit
        status.
    */
    static bool runCmd(const std::string& cmdLine,
                       const int expectedExitStatus = 0);
    
protected:
    // Currently this class does not have a protected member.
    
private:
    /** Default constructor.

        This class is not meant to be instantiated. Instead the
        various static helper methods can be directly called to
        perform various operations.
    */
    Helpers() {}

    /** Destructor

        The destructor is private to ensure that an object of this
        type can never be deleted (at it cannot be created either).
    */
    ~Helpers() {}
};

#endif
