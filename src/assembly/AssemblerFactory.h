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
class ArgParser;

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
    friend class AssemblySubSystem;
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

        \param[in] mpiRank The MPI assigned rank for the process on
        which the assembler is being instantiated.  Some of the
        assemblers use this information to further customize their
        operations.
    */
    static Assembler* create(const std::string& name,
                             const int mpiRank);

protected:
    /** A set of ArgParser::ArgRecord entries to display supported
        heuristics.

        This static list of heuristic names is used by the
        ClusteringSubSystem::addCommandLineArguments() method to
        conveniently add the list the supported set of valid
        analyzers.  The functionality of listing the valid analyzers
        is present here because this class is responsible for
        instantiating ESTAnalyzer classes and its the one that is
        aware of the various analyzers.
		
        \param[out] argParser The argument parser to which the list of
        valid entries are to be added in the form of
        ArgParser::INFO_MESSAGE.
    */
    static void addCommandLineInfo(ArgParser& argParser);
    
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
