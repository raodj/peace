#ifndef ASSEMBLER_FACTORY_H
#define ASSEMBLER_FACTORY_H

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

#include <iostream>

/** \file AssemblerFactory.h
    
    \brief A convenient single-point object creation factory for all
    gene assemblers in PEACE.

    This file contains the class declaration for the
    AssemblerFactory. The AssemblerFactory must be used to instantiate
    a specific gene assembler from the family of related assembler
    classes.  Additional documentation regarding the API provided by
    this class is available with the class documentation.
*/

// Forward declaration to make compiler happy and keep it fast
class Assembler;

/** An object factory for instantiating a specific assembler.

    This class serves as a single-point interface for instantiating
    one of the different types of assemblers supported by PEACE.  This
    class has been implementated in concordance with the standard
    object-oriented computer programming strategy of an object
    factory.  This class is primarily used to create derived assembler
    classes that are child classes of the Assembler API class.  Note
    that this class is a non-instantiable class.  However, its API
    methods are static and consequently they can be directly invoked
    without requiring an instance.
*/
class AssemblerFactory {
public:
    /** Method to instantiate a suitable Assembler object.

        This method must be used to instantiate a suitable assembler
        instance.  This method uses the name (parameter) to suitably
        instantiate an assembler.  If the name is not valid, then this
        method returns NULL.

        \note The Message Passing Interface (MPI), if being used, must
        have been initialized before this method is invoked.
        
		\param[in] name The name of the derived gene assembler to be
        instantiated.

        \param[in] outputFileName The target file to which the
        assembled gene sequences and other information is to be
        written (if any).  If this parameter is the empty string (""),
        then the output is written to standard output.

        \param[in] mpiRank The MPI assigned rank for the process on
        which the assembler is being instantiated.  Some of the
        assemblers use this information to further customize their
        operations.
    */
    static Assembler* create(const char* name,
                             const std::string& outputFileName,
                             const int mpiRank);

    /** Method to display the list of assemblers available.

        This method is typically used to display the list of
        assemblers currently compiled into the PEACE system.  This
        method is typically used in the main() method when displaying
        usage information.

        \param[out] os The output stream to which the list of EST
        names must be written.
    */
    static void displayList(std::ostream& os);
    
protected:
    // Currently this class has no protected members.
    
private:
    /** The default constructor.

        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        in this class must be directly used to instantiate a cluster
        maker.
    */
    AssemblerFactory() {}

    /** The destructor.

        The destructor is private to ensure that objects of this class
        are never deleted.
    */
    ~AssemblerFactory() {}
};


#endif
