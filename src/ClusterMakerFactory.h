#ifndef CLUSTER_MAKER_FACTORY_H
#define CLUSTER_MAKER_FACTORY_H

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
class ClusterMaker;
class ESTAnalyzer;

class ClusterMakerFactory {
public:
    /** Method to instantiate a suitable Cluster Maker object.

        This method must be used to instantiate a suitable cluster
        maker instance.  This method uses the name (parameter) to
        suitably instantiate an cluster maker.  If the name is not
        valid, then this method returns NULL.

	\param[in] name The name of the cluster maker to be
        instantiated.
	
	\param[inout] analyzer The EST analyzer to be used by this
        ClusterMaker for generating similarity metrics between two
        given ESTs. If this parameter is NULL, this method prints an
        error message and returns NULL.
	
        \param[in] refESTidx The reference EST index value to be used
        when instantiating a cluster maker.  Index values start from 0
        (zero).  If this value is negative then this method returns
        NULL.

        \param[in] outputFileName The target file to which the
        analysis report is to be written (if any).  If this parameter
        is the empty string (""), then the output is written to
        standard output.
    */
  static ClusterMaker* create(const char* name, ESTAnalyzer *analyzer,
				const int refESTidx,
                                const std::string& outputFileName);

    /** Method to display the list of cluster makers available.

        This method is typically used to display the list of cluster
        makers currently compiled into the program.  This method is
        typically used in the main() method when displaying usage
        information.

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
    ClusterMakerFactory() {}

    /** The destructor.

        The destructor is private to ensure that objects of this class
        are never deleted.
    */
    ~ClusterMakerFactory() {}
};


#endif
