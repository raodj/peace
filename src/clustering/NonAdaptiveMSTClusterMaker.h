#ifndef NON_ADAPTIVE_MST_CLUSTER_MAKER_H
#define NON_ADAPTIVE_MST_CLUSTER_MAKER_H

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

#include "MSTClusterMaker.h"

/** A Minimum Spanning Tree (MST) based parallel cluster maker.

    This class encapsulates the core functionality needed to construct
    a MST-based EST clusters in a parallel/distributed manner using
    the Message Passing Interface (MPI) library.  This class includes
    functionality for both the Manager (MPI Rank == 0) and Worker (MPI
    Rank > 0) processes.  Necessary functionality to distinguish and
    operate either as Manager or Worker is already built into the
    class.  This class uses the MSTCache and MSTCluster classes to
    help in performing the various activities.  Refer to the
    documentation on the various method for detailed description on
    their functionality and usage.
*/
class NonAdaptiveMSTClusterMaker : public MSTClusterMaker {
    friend class ClusterMakerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~NonAdaptiveMSTClusterMaker();

protected:
    // Currently this class does not have protected members

private:
    /** The constructor.
        
        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ClusterMakerFactory is a
        friend of this class, an object can be instantiated via the
        ClusterMakerFactory::create() method.

        \param[in,out] analyzer The EST analyzer to be used for
        obtaining similarity metrics between two ESTs.  This parameter
        is simply passed onto the base class.
    */
    NonAdaptiveMSTClusterMaker(ESTAnalyzer *analyzer);    
};

#endif
