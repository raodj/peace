#ifndef BIPARTITE_DATA_H
#define BIPARTITE_DATA_H

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

/** A class to encapsulate data and provide convenience methods for
    obtaining information on a bipartite graph.

    Bipartites are used to connect the separate partitions in the
    PMSTClusterMaker algorithm.  A bipartite consists of the ESTs from
    two partitions, with the catch that ESTs from one partition will
    only be compared with ESTs from the other partition in constructing
    the MST.  Additionally complicating matters is the fact that, since
    bipartite MSTs take longer to construct than MSTs on single partitions,
    multiple processes will be assigned to work on a bipartite whenever
    possible.  This class attempts to encapsulate information on the
    bipartite as well as the processes working on it and provide methods
    for those processes to conveniently access relevant data (such as
    the rank of the process acting as manager for building the MST).

    \note This class is meant to be used only by the PMSTClusterMaker.
*/
class BipartiteData : public PartitionData {
public:
    /** The constructor.
      
        This is a convenience constructor that is used to set up all
        the instance variables to corresponding values specified by
        the parameters.

        \param[in] startESTidx The zero-based index of the EST that
	marks the beginning of this partition.
          
        \param[in] estCount The count of ESTs belonging to this partition.
    */
    BipartiteData(const int startESTidxValue, const int estCountValue,
		  PartitionData* p1, PartitionData* p2,
		  const int managerRank, const int workerCount) :
      PartitionData(startESTidxValue, estCountValue,
		    estCountValue/(workerCount+1)), partition1(p1),
      partition2(p2), manager(managerRank), workers(workerCount) {}
    
    /** The destructor.

        This class does not contain any members that utilize dynamic
        memory to hold data.  Consequently, the destructor has no
        special tasks to perform and is merely present to adhere to
        coding conventions.
    */
    virtual ~BipartiteData() {
      delete partition1;
      delete partition2;
    }

    //inline void getOwnedESTs(const int rank, int* startIdx1, int* startIdx2) {
    //  
    // }

    virtual int getPartitionManager() const
      { return manager; }

    virtual int getWorkerCount() const
      { return workers; }

    /** A pointer to the first partition in this bipartite.

    */
    PartitionData* partition1;

    /** A pointer to the second partition in this bipartite.

    */
    PartitionData* partition2;

    /** The MPI rank of the manager for this bipartite.

    */
    int manager;
    
    /** The count of worker processes assigned to this bipartite.

    */
    int workers;

protected:
    // Currently this class does not have any protected members.
    
private:
    // Currently this class does not have any private members.
};

#endif
