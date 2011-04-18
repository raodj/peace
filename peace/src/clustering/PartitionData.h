#ifndef PARTITION_DATA_H
#define PARTITION_DATA_H

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

/** A simple class to encapsulate data regarding a partition in the
    PMSTClusterMaker owned by a particular process.  A partition is a
    subset of the EST data set.

    \note This class is meant to be used only by the PMSTClusterMaker.
*/
class PartitionData {
public:
    /** The constructor.

        This constructor simply sets all the instance variables to
		corresponding values specified by the parameters.
		
        \param[in] startESTidxValue The zero-based index of the EST
		that marks the beginning of this partition.
		
        \param[in] estCountValue The count of ESTs belonging to this
        partition.
    */
    PartitionData(const int startESTidxValue, const int estCountValue) :
        startESTidx(startESTidxValue), estCount(estCountValue),
        ownedESTCount(estCountValue) {}
    
    /** The constructor.

        This constructor simply sets all the instance variables to
		corresponding values specified by the parameters.

        \param[in] startESTidxValue The zero-based index of the EST
        that marks the beginning of this partition.
        
        \param[in] estCountValue The count of ESTs belonging to this
        partition.
        
        \param[in] ownedESTCountValue The count of ESTs owned by this
        partition.
    */
    PartitionData(const int startESTidxValue, const int estCountValue,
                  const int ownedESTCountValue) :
        startESTidx(startESTidxValue), estCount(estCountValue),
        ownedESTCount(ownedESTCountValue) {}
    
    /** The destructor.

        This class does not contain any members that utilize dynamic
        memory to hold data.  Consequently, the destructor has no
        special tasks to perform and is merely present to adhere to
        coding conventions.
    */
    virtual ~PartitionData() {}

    virtual int getPartitionManager() const
      { return -1; }

    virtual int getWorkerCount() const
      { return 0; }

    /** The zero-based index of the first EST in the partition.

    */
    int startESTidx;
    
    /** The count of ESTs belonging to the partition.

    */
    int estCount;

    /** The count of ESTs which are owned by this process.

    */
    int ownedESTCount;

protected:
    // Currently this class does not have any protected members.
    
private:
    // Currently this class does not have any private members.
};

#endif
