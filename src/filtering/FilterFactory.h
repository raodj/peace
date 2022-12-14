#ifndef FILTER_FACTORY_H
#define FILTER_FACTORY_H

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

#include <string>

// Forward declarations to make compiler happy and keep it fast.
class Filter;
class RuntimeContext;

/** Centralized factory for creating filters.

	This is a centralized helper class that can be used to create new
	instances of various filters. This class provides a convenient
	FilterFactory::create() method for this purpose. In addition, it
	also provides a mechanism to list all the known filters.

    <p>A set of filters (stored in the FilterChain class) are run on
    all entries in RuntimeContext::estList prior to commencement of
    the core clustering/assembly operation. The filters perform
    various validation operations to ensure ESTs are good prior to
    clustering. Filtering ensures that the overall quality of
    clustering provided by PEACE is good.</p>

    <p>Each filter object in the chain implements a specific type of
    filteration operation and ultimately returns an integer indicating
    the cluster to which an EST is to be assigned. If the cluster ID
    is -1, then that indicates that the EST must be subjected to
    regular clustering operations.</p>

    \see FilterChain
*/
class FilterFactory {
public:
    /** Method to instantiate a suitable filter.
        
        This method must be used to instantiate a suitable filter
        instance.  This method uses the name (parameter) to suitably
        instantiate a filter.  If the name is not valid, then this method
        returns NULL.
        
        \param[in] name The name of the filter to be instantiated.

        \param[in] runtimeContext The runtime context class from which
        the necessary context information (such as ESTAnalyzer to use
        or ClusterMaker to use) is to be extracted.  The pointer is
        saved so that the filter class can use it to perform specific
        operations.
    */
    static Filter* create(const std::string& name,
                          RuntimeContext *runtimeContext);   
protected:
    // Currently this class has no protected members.
    
private:
    /** The default constructor.
        
        The default constructor has is private in order to ensure that
        this class is never instantiated.  Instead, the static methods
        in this class must be directly used to instantiate a filter.
    */
    FilterFactory() {}
    
    /** The destructor.
        
        The destructor is private to ensure that objects of this class
        are never deleted.
    */
    ~FilterFactory() {}
};

#endif
