#ifndef MST_CACHE_H
#define MST_CACHE_H

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

#include <vector>
#include <utility>
#include <functional>
#include <iostream>

#include "ESTAnalyzer.h"
#include "CachedESTInfoHelper.h"

/** \def SMList

    \brief Typedef for std::vector<CachedESTInfo>.

    This typedef is a convenience definition to handle a vector of
    CachedESTInfo objects associated with a given EST in the MST
    Cache.  This data type is used by the MSTClusterMaker to create a
    temporary list of cache entries and then add a whole bunch of
    entries to the MSTCache.  Adding a bunch of cache entries at a
    time facilitates distribution of data (to various workers) and
    streamlines management of cached information.
*/
typedef std::vector<CachedESTInfo> SMList;

/** A cache to maintain similarity/distance metrics to build a Minimum
    Spanning Tree (MST).

    <p>This class is the base class of all distributed-cache
    maintainers used for rapdily building a Minimum Spanning Tree
    (MST) for clustering.  This class defines the core API for all
    ESTCache classes.  In addition, it includes a few helper methods
    that provides some convenience methods to manage the caches.
    Furthermore, this class streamlines the operation of the Master
    and Worker processes when building a MST in a distributed manner.
    This class (and its children) essentially cache
    similarity/distance metrics to minimize the number of times the
    expensive operation of computing similarity metrics is performed
    via an EST analyzer. </p>
*/
class MSTCache {
public:
    /** The destructor.
        
        The parent class (std::vector) manages all the memory
        operations for this class.  Consequently, the destructor does
        not have any special tasks to perform.  It here present merely
        to adhere to coding conventions.
    */
    virtual ~MSTCache() {}
    
    /** Determine if a EST (given its index) is already present in the
        MST.

        This method can be used to determine if a given EST has been
        added to the MST.  This method was introduces to make the code
        a bit more readable at certain spots.

        \param[in] estIdx The index of the EST to be tested.

        \return This method returns \c true if the EST (with the given
        index) is present in the MST.  Otherwise it returns \c false.
    */
    inline bool isESTinMST(const int estIdx) const
    { return nodesAdded[estIdx]; }

    /** Suggest to clear caches that contain entries for a newly added
        EST.

        This method can be used to suggest to the cache to prune
        entries corresponding to the \c estIdxAdded from the various
        data structures maintained by this cache.  The pruning
        operation is useful for the following two purposes:

        <ol>

        <li>First it removes entries from various lists corresponding
        to estIdxAdded (because these entries are vestigial as the
        node as already been added to the MST) thereby reducing memory
        usage. </li>

        <li>Second, if a list becomes empty then the cache needs to be
        repopulated with fresh entries for further MST
        computations.</li>

        </ol>

        \note Some dervied classes may choose not to prune caches or
        add entries to the \c repopulateList. Consequently, this
        method must be viewed as a suggestion to the \c MSTCache to
        give it a chance to prune caches as needed.
        
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
                             const bool prefixListSize = true) = 0;

    /** Give the MSTCache a chance to preprocess lists before
        distributing them over the network.

        This method is a helper method that can be used by the
        MSTClusterMaker to have a local metric list to be preprocessed
        prior to dispatching it over the network. Preprocessing a
        local list possibly reduces the size and distributes some of
        the processing overheads (thereby accelerating overall
        processing).

        \note The default method implementation in the base class does
        absolutely nothing.
        
        \param[in,out] list The list of entries that must be
        preprocessed.
    */
    virtual void preprocess(SMList& UNREFERENCED_PARAMETER(list)) {}
    
    /** Add/merges a batch of new cache entries with the current
        entries in the cache for a given EST.

        This method merges a given batch of new cache entries with
        existing cache entries for the given EST.  Some dervied
        classes may ensure that only MaxCacheSize entries are retained
        in the cache (retained entries \b must be chosen so that they
        have the best metrics).  Entries in the cache are removed to
        accommodate new entries.

        \note Derived classes must override this method.  However,
        they may choose to ignore the \c MaxCacheSize suggestion and
        retain more entries than necessary to reduce unnecessary
        computations.
        
        \param[in] estIdx The zero-based index of the reference EST
        against which the batch of cache entries in \c list were
        analyzed and generated.

        \param[in] list A SMList containing the batch of new cache
        entries to be merged with any existing entries in the cache
        for the given EST (whose zero-based index is \c estIdx).
    */
    virtual void mergeList(const int estIdx, const SMList& list) = 0;

    /** Obtains the top-most similar entry from the MSTCache.

        This method searches the similarity metrics stored in this
        cache for various ESTs to locate the best entry with the
        highest similarity.  It populates the parameters with the
        appropriate value.  Note that the parameters are initialized
        to -1, -1, -1.0f, and -1 respectively.

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
			      int& directionData) const = 0;

    /** Display cache usage statistics.

        This method can be used to print the current cache usage
        statistics.  Different derived classes track and report
        different statistics.  Consequently, the actual statistics
        reported by this method will vary.

        \note Derived classes must override this method.
        
        \param[out] os The output stream to which statistics must be
        written.
        
        \param[in] MyRank The MPI rank of the process to which this
        cache is assocaited.  This information is used merely for
        displaying more informative statistics and does not influence
        the actual values displayed by this method.
    */
    virtual void displayStats(std::ostream &os, const int MyRank) const = 0;
    
protected:
     /** The constructor.

        The MST cache requires information about the total number of
        ESTs in the system in order to effectively track the ESTs
        already added to the MST.  This value must be passed in as the
        only parameter.

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
    MSTCache(const int totalESTCount, const int startOwnESTidx,
             const int numOwnESTs, const ESTAnalyzer *analyzer,
             const bool repopulateCaches, const int maxCachePerEST);

    /** Utility method to copy first n-entries from input to output SMList.

        This is a rather straightforward method that can be used to
        copy the first \c cout entries from the \c input list and
        append them to the \c output list.

        \param[in] input The list of entries that must be copied to
        the \c output list.

        \param[in] count The number of entries to be copied.  If this
        value is greater then input.size(), then only input.size()
        entries are copied.

        \param[out] output The output list to which the entries are to
        be copied.
    */
    static void copy_n(const SMList& input, const size_t count, 
                       SMList &output);

    /**  The analyzer to be used for any EST comparison.
         
         This pointer is initialized in the constructor to point to
         the EST analyzer that is used for EST comparisons.  This
         pointer is essentially used for comparing EST metrics for
         ordering information in the cache.
    */
    const ESTAnalyzer* const analyzer;
        
    /** Bit-vector to determine if a given node has been added to the
        MST.
        
        This vector is used to track if a given EST node has been been
        added to the MST.  This std::vector<bool> is a specialization
        that optimizes memory usage.  So this member should not occupy
        much memory.  The maximum number of nodes is set in the
        constructor and all entries are initialized to false.  The
        pruneCache() method sets the suitable entry to \c true after
        the EST has been added to the MST.  The isESTinMST() method
        uses this vector.
    */
    std::vector<bool> nodesAdded;

    /** The index of the first EST whose information is cached.

        The starting index of the EST that is owned by the process
        that is using this cache.  This value is used by some derived
        classes to internally normalize EST index values to zero for
        the first EST owned by this process.  This value is
        initialized in the constructor and is never changed during the
        lifetime of this object.
    */
    const int startOwnESTidx;

    /** The number of ESTs that this cache logically owns and manages.
        
        This instance variable maintains the number of ESTs that are
        logically owned and managed by this cache.  This value is
        initialized in the constructor and is never changed during the
        lifetime of this object.
    */
    const int numOwnESTs;

    /** The \b suggested maximum number of cached entries per EST.

        This instance variable is a \b suggested maximum number of
        cached entries that must be maintained for a given EST. This
        parameter may not be really used by all dervied classes (as
        this is only a suggestion).
    */
    const int maxCachePerEST;
    
    /** Flag to control if caches must be repopulated.

	<p> If this parameter is \c true then the MSTCache will
	request lists to be repopulated when needed. Repopulating
	lists guarantees that ultimately a MST will be developed.  If
	repopulation is suppressed then the resulting spanning tree
	may not be a MST. </p>

        <p> This value is set in the constructor and is never changed
        during the life time of this object. </p>
        
    */
    const bool repopulateCaches;

    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    MSTCache& operator=(const MSTCache& src);
};

#endif
