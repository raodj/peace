#ifndef HEURISTIC_FACTORY_H
#define HEURISTIC_FACTORY_H

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

// For declaration to make compiler happy.
class Heuristic;
class HeuristicChain;
class ArgParser;

/** Centralized factory for creating heuristics.

    This is a centralized helper class that can be used to create new
    instances of various heuristics. This class provides a convenient
    HeuristicFactory::create() method for this purpose.
    
    <p>A set of heuristics (stored in the HeuristicChain class) are
    run on selected entries in RuntimeContext::estList prior to
    running heavy weight clustering/assembly operation. The heuristics
    perform various validation operations to ensure that only
    sufficiently similar ESTs are subject to detailed
    clustering/assembly.  The heuristics run much faster than the core
    heavy-weight algorithms and thereby improve performance.</p>

    <p>Each heuristic object in the chain implements a specific type
    of heuristic operation and ultimately returns a boolean value
    indicating if two fragments are sufficiently similar to waranning
    running the heavy weight algorithms on them.</p>

    \see HeuristicChain
*/
class HeuristicFactory {
    friend class ClusteringSubSystem;
public:
    /** Method to instantiate a suitable heuristic.
        
        This method must be used to instantiate a suitable heuristic
		instance.  This method uses the name (parameter) to
        suitably instantiate a heuristic.  If the name is not
        valid, then this method returns NULL.
        
        \param[in] name The name of the heuristic to be instantiated.

        \param[in] chain The heuristic chain is going to logically
        contain this heuristic.  This value is used by the heuristic
        object to update hints (used by heavy weight analyzers) in the
        heuristic chain class.  This value is saved by the heuristic
        object and is never changed or deleted by the Heuristic class
        (or its children).
    */
    static Heuristic* create(const std::string& name, HeuristicChain *chain);
    
protected:
    /** A set of ArgParser::ArgRecord entries to display supported
        heuristics.

        This static list of heuristic names is used by the
        ClusteringSubSystem::addCommandLineArguments() method to
        conveniently add the list the supported set of valid
        heuristics.  The functionality of listing the valid heuristics
        is present here because the HeuristicFactory is responsible
        for instantiating heuristics and its the one that is aware of
        the various heuristics.

        \param[out] argParser The argument parser to which the list of
        valid entries are to be added in the form of
        ArgParser::INFO_MESSAGE.
    */
    static void addCommandLineInfo(ArgParser& argParser);
    
private:
    /** The default constructor.
        
        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        in this class must be directly used to instantiate a heuristic.
    */
    HeuristicFactory() {}

    /** The destructor.

        The destructor is private to ensure that objects of this class
        are never deleted.
    */
    ~HeuristicFactory() {}
};

#endif
