#ifndef MST_MULTI_LIST_CACHE_H
#define MST_MULTI_LIST_CACHE_H

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

/** \typedef std::pair<int, SMList> MSTCacheEntry

    \brief Shortcut typedef for std::pair<int, SMList>

    This typedef is a shortcut for a cache entry that contains the
    following pair of information:

    <ol>

    <li>The first entry is an \c int that indicates the EST index with
    which this cache entry is associated. </li>

    <li>The second entry is a \c SMList that contains the closet (most
    similar) ESTs for the EST indicated by the first entry.</li>

    </ol>
*/
typedef std::pair<int, SMList> MSTCacheEntry;

/** A cache that uses multiple lists to maintain similarity/distance
    metrics to build a Minimum Spanning Tree (MST).

    This class dervies off the \c MSTCache base class provides the
    functionality to manage a cache of similarity metrics using a
    collection of Lists.  This cache is organized as follows: The
    cache is essentially a collection of SMList objects, one SMList
    for each EST (logically) owned by this cache.  The lists are
    maintained in a sorted fashion with the best metrics at the top of
    the list to facilitate rapid construction of a Minimum Spanning
    Tree (MST).
*/
class MSTMultiListCache : public MSTCache {
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
    MSTMultiListCache(const int totalESTCount, const int startOwnESTidx,
                      const int numOwnESTs, const ESTAnalyzer *analyzer,
                      const bool repopulateCaches, const int maxCachePerEST);

    /** The destructor.
        
        The instance variables in this class automatically manage all
        the memory associated this class.  Consequently, the
        destructor does not have any special tasks to perform.  It
        here present merely to adhere to coding conventions.
    */
    virtual ~MSTMultiListCache() {}

    /** Clear caches that contain entries for a newly added estIdx.

        This method overrides the virtual method in the base class and
        prunes entries corresponding to the \c estIdxAdded from the
        various lists maintained by this cache.  The pruning operation
        is removes unused entries and facilitates repopulation of
        caches (if requested).

        \param[in] estIdxAdded The index of the EST that has just been
        added to the MST.

        \param[out] repopulateList The index values of ESTs whose
        cache needs to be repopulated.  The list of EST indexes added
        to this list are only the ESTs that are owned by this process.
        If none of the caches need to be recomputed then this list is
        empty (all existing entries are removed).

	\param[in] prefixListSize If this parameter is set to \c true,
	then the number of entries added to the repopulateList is
	inserted at the beginning of the vector.  This eases
	transmission of lists via MPI to the master process.
    */
    virtual void pruneCaches(const int estIdxAdded,
                             std::vector<int>& repopulateList,
                             const bool prefixListSize = true);

    /** Merges a given SMList with current entries in the cache for a
        given EST.

        This method overrides the virtual method in the base class. As
        per the API requirement, it merges the given list with
        existing cache entries for the given EST such that only \c
        maxCachePerEST entries are retained in the cache.  The
        retained entries are chosen so that they have the best
        metrics.  Care is taken to ensure that the cache continues to
        remain sorted with the best metrics at the top of the cache.

        \param[in] estIdx The index of the EST with which the given
        list must be merged.

        \param[in] list The SMList to be merged with any existing
        entries in the cache for the given estIdx.
    */
    virtual void mergeList(const int estIdx, const SMList& list);

    /** Obtains the top-most similar entry from the MSTCache.

        This method overrides the default implementation in the base
        class.  It searches the analysis metrics stored in this cache
        for various ESTs to locate the best entry with the best
        metric.  It populates the parameters with the appropriate
        value.  Note that the parameters are initialized to -1, -1,
        and -1.0f respectively.

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
    */
    void getBestEntry(int& srcESTidx, int& destESTidx,
                      float& metric, int& alignmentData) const;


    /** Display cache usage statistics.

        This method can be used to print the current cache usage
        statistics.  This method prints the following information:

        <ul>

        <li>Number of ESTs cached by this cache.</li>
        
        <li>Total number of SMEntries currently in the cache.</li>

        <li>Average length of each list and the standard deviation.</li>

        <li>The number of times the cache had to be repopulated.</li>

        </ul>

        \param[out] os The output stream to which statistics must be
        written.

        \param[in] MyRank The MPI rank of the process to which this
        cache is assocaited.  This information is used merely for
        displaying more informative statistics and does not influence
        the actual values displayed by this method.
    */
    void displayStats(std::ostream &os, const int MyRank) const;
    
    /** Sort a Similarity Metric List (SMList) based on similarity
        metric and then prunes it.
        
        This method provides a convenient method for sorting a
        similarity metric list based on the similarity metric.  This
        method simply uses the STL sort method with a suitable Functor
        for performing this task.  Sorting is performed such that the
        highest similarity metric is first in the vector after
        sorting.  Once sorting has been completed, this method ensures
        that the list is no longer than the specified cacheSize.
 
        \param[inout] list The SMList to be sorted based on the
        (similarity or distance) metric associated with each
        CachedESTInfo entry in the list.

	\param[in] cacheSize The maximum size of this list after it
	has been sorted.
    */
    void preprocess(SMList& list);

protected:    
    /** Helper method to prune a given SMList.

        This is a helper method that is called from the pruneCaches()
        method to remove all entries corresponding to the given
        estIdx.

        \param[inout] list The SMList whose entries needs to be
        pruned.

        \param[in] estIdx The index of the EST whose entries must be
        removed from this list.

        \return This method returns \c true if the list becomes empty
        once all the entries have been removed.
    */
    bool pruneCache(SMList& list, const int estIdx);

private:   
    /** The actual vector of cache entries for various nodes.

        This cache contains the list of cache entries for each node
        owned and managed by this MSTCache.  The list is initialized
        to have a set of empty lists in the constructor, when the
        cache is instantiated for use.
    */
    std::vector<MSTCacheEntry> cache;

    /** Number of times the cache needed to be repopulated.
 
        This instance variable is used to track the number of times
        the cache had to be repopulated as one of the entries ran out
        of data.  This variable is initialized to 0 (zero),
        incremented in the pruneCaches() method, and displayed in the
        displayStats() method.
    */
    int cacheRepopulation;

    /** Number of entries that were pruned from the cache.

        This instance variable is used to track the number of entries
        that were pruned from the cache by the pruneCache() method.
        This value indicates the number of entries in the cache that
        were unused.  Moreover, this value also a measure of the
        amount of unnecessary adjacency information calculation
        operations that were performed for this cache by all the
        processes put together.
    */
    int prunedEntries;

    /** Number of times comparisons of ESTs were performed.
        
        Tracks the number of times EST comparisons were performed.
        This instance variable is initialized to 0 (zero).  The set of
        values sorted via teh sortAndPrune() method are used to track
        the number of comparisons performed and this instance variable
        is suitably adjusted.
    */
    long analysisCount;    
};

#endif
