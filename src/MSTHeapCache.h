#ifndef MST_HEAP_CACHE_H
#define MST_HEAP_CACHE_H

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

#include "MSTCache.h"
#include <queue>

/** A cache to maintain similarity/distance metrics to build a Minimum
    Spanning Tree (MST) using a heap.

    This class dervies off the \c MSTCache base class provides the
    functionality to manage a cache of metrics using a single large
    heap.  The operations on a heap are fast. However, the current
    heap does not allow to rapidly prune unwanted entries right
    away. Consequently, pruning of entries is done on a
    "whenever-chance-arises" basis. This causes \c MSTHeapCache to
    have a much higher memory foot print when compared to
    MSTMultiListCache.
*/
class MSTHeapCache : public MSTCache {
public:
    /** The constructor.
        
        This MST cache requires information about the number of ESTs
        in the system in order to effectively track the ESTs already
        added to the MST.  These values must be passed in as the
        various parameters to this class.  The constructor essentially
        passes off the parameters to the base class that saves the
        information (for future reference) in suitable instance
        variables.

        \note This object is typically instantiated in the
        MSTClusterMaker::makeClusters() method.

        \param[in] totalESTCount The total number of ESTs to be
        processed (that is, to be added to the MST).

        \param[in] startOwnESTidx The starting index of the EST that
        is owned by the process that is using this cache.  This value
        is used to internally normalize est index values to 0 for the
        first EST owned by this process.

        \param[in] numOwnESTs The number of ESTs owned by the process
        that is using this cache.  This information is used to reserve
        necessary space to hold SMLists for each EST owned by the
        manager/worker process.

        \param[in] analyzer The analyzer to be used for EST
        comparisons, comparing metrics, and ordering entries in this
        cache depending on the type of EST analyzer (whether it is
        similarity based or distance metric based).

	\param[in] repopulateCaches If this parameter is \c true then
	the MSTCache will request lists to be repopulated when
	needed. Repopulating lists guarantees that ultimately a MST
	will be developed.  If repopulation is suppressed then the
	resulting spanning tree may not be a MST.

        \param[in] maxCachePerEST This parameter \b suggests the
        maximum number of cached entries that must be maintained for a
        given EST. This parameter may not be really used by all
        dervied classes (this is only a suggestion).
    */
    MSTHeapCache(const int totalESTCount, const int startOwnESTidx,
                 const int numOwnESTs, const ESTAnalyzer *analyzer,
                 const bool repopulateCaches, const int maxCachePerEST);
    
    /** The destructor.
        
        The STL data structures used by this class manage all the
        memory operations for this class.  Consequently, the
        destructor does not have any special tasks to perform.  It
        here present merely to adhere to coding conventions.
    */
    virtual ~MSTHeapCache() {}
    
    /** Suggest to clear caches that contain entries for a newly added
        EST.

        This method can be used to suggest to the cache to prune
        entries corresponding to the \c estIdxAdded from the various
        data structures maintained by this cache. This method
        repeatedly removes entries at the top of the heap if the
        entries correspond to estIdx.
        
        \param[in] estIdxAdded The index of the EST that has just been
        added to the MST.

        \param[out] repopulateList The index values of ESTs whose
        cache needs to be repopulated.  Currently, this method does
        not populate any entries into this list.

	\param[in] prefixListSize If this parameter is set to \c true,
	then the number of entries added to the repopulateList is
	inserted at the beginning of the vector.  This eases
	transmission of lists via MPI to the master process.
    */
    virtual void pruneCaches(const int estIdxAdded,
                             std::vector<int>& repopulateList,
                             const bool prefixListSize = true);
   
    /** Add/merges a batch of new cache entries with the current
        entries in the cache for a given EST.

        This method merges a given batch of new cache entries with
        existing cache entries in the cache.  Caches are not trimmed
        to fit into the suggested maxCachePerEST size.
        
        \param[in] estIdx The zero-based index of the reference EST
        against which the batch of cache entries in \c list were
        analyzed and generated.

        \param[in] list A SMList containing the batch of new cache
        entries to be merged with any existing entries in the cache
        for the given EST (whose zero-based index is \c estIdx).
    */
    virtual void mergeList(const int estIdx, const SMList& list);

    /** Obtains the top-most similar entry from the MSTCache.

        This method simply uses the top-most entry in the heap and
        populates the parameters with the appropriate value.  Note
        that the parameters are initialized to -1, -1, -1.0f, and -1
        respectively.

        \param[out] srcESTidx The source EST index from where the
        similarity metric is being measured.  The srcESTidx is already
        present in the MST.

        \param[out] destESTidx The destination EST index that is the
        best choice to be added to the MST (based on the local
        information).

        \param[out] metric The similarity/distance metric between the
        srcESTidx and the destESTidx.

        \param[out] alignmentData The alignment data associated with
        the srcESTidx and the destESTidx.

	\param[out] directionData The direction data associated with
	the srcESTidx and the destESTidx.
    */
    virtual void getBestEntry(int& srcESTidx, int& destESTidx,
                              float& metric, int& alignmentData,
			      int& directionData) const;

    /** Display cache usage statistics.

        This method currently displays the following information:

        <ul>
        
        <li>Total number of entries cached by this cache.</li>
        
        <li>Current number of entries left in the cache.</li>
        
        </ul>
        
        \param[out] os The output stream to which statistics must be
        written.
        
        \param[in] MyRank The MPI rank of the process to which this
        cache is assocaited.  This information is used merely for
        displaying more informative statistics and does not influence
        the actual values displayed by this method.
    */
    virtual void displayStats(std::ostream &os, const int MyRank) const;
    
protected:
    /** A simple wrapper method to handle popping the cache (removing
	the top element from the cache) and incrementing the number
	of pruned entries.
    */	
    virtual void popCache();
    
private:
    /** Number of entries that were pruned from the cache.

        This instance variable is used to track the number of entries
        that were removed from the cache by the pruneCache() method.
    */
    int prunedEntries;

    /** The actual heap of cache entries for various nodes.

        This cache contains all the cached entries for each EST owned
        and managed by this MSTCache. Entries are added by the \c
        mergeList method and removed by the \c pruneCache method.
    */    
    std::priority_queue<CachedESTInfo, std::vector<CachedESTInfo>,
                        GreaterCachedESTInfo> cache;

    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    MSTHeapCache& operator=(const MSTHeapCache& src);
};

#endif
