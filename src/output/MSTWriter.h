#ifndef MST_WRITER_H
#define MST_WRITER_H

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

#include "Component.h"

// Forward declarations to keep compiler fast and happy
class MST;

/** An output sub-system component to write Minimum Spanning Tree
	(MST) data.

	This class is a component of the output sub-system. Its
	responsibility is to serialize a given Minimum Spanning Tree (MST)
	data structure.  This class permits serialization of MST data
	either to a given data file.

	\note This class can only be instantiated by the
	OutputSubSystem. See OutputSubSystem::getMSTWriter() method.
*/
class MSTWriter : public Component {
    friend class OutputSubSystem;
public:
    /** The destructor.

        The destructor is merely present to adhere to coding
        conventions.  Currently, this does not have any special tasks
        to perform.
    */
    virtual ~MSTWriter();

    /** Add the set of command line parameters for this component.

        This method is typically invoked by the OutputSubSystem when
        it is being requested for command line arguments.  This method
        adds the \c --output-mst-file option to the set of command
        line arguments that can be used to further customize the
        operations of this component.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Write the given MST data.

        This method must be used to appropriately serialize the MST
        data as directed by the user (through command-line arguments).

        \param[in] mst The Minimum Spanning Tree (MST) to be
        serialized by this method.  

        \return If the MST is written successfully, then this method
        returns \c true.  On errors it generates a suitable error
        message (on std::cerr) and returns \c false.
    */
    virtual bool write(const MST& mst);
    
protected: 
    /** The constructor.

        This is the only constructor for this class. The
        constructor(s) are not public to ensure that this class is not
        instantiated directly but only by the friends of this
        component.  Currently, only the OutputSubSystem class can
        instantiate this component.  The constructor does not have any
        special tasks to perform and merely initializes instance
        variables to their default initial value.
    */
    MSTWriter();

private:
    /** The name of the file to which the MST data is to be written.

        This instance variable contains the path (may it be relative
        or absolute) and name of the file to which the MST data is to
        be written when the write() method in this class is invoked.
    */
    std::string mstFileName;
};

#endif
