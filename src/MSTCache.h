#ifndef MST_CACHE_H
#define MST_CACHE_H

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include <vector>
#include <utility>
#include <functional>
#include <iostream>

#include "ESTAnalyzer.h"

/** \typedef SMEntry

    \brief Typedef for an entry in the Similarity/Distance Metric (SM)
    List (std::pair<int, float>).
    
    A Similarity/Distance Metric (SM) Entry is associated with a given
    EST.  The EST with which this SMEntry is associated is \b not
    stored in the SMEntry to save meory.  Instead, a
    std::vector<SMEntry> (aka SMList) is explicitly associated with a
    given EST index.  An SMEntry consists of the following two pieces
    of information:

    <ol>

    <li> An \c int representing the EST index. </li>

    <li> A \c float representing the similarity/distance metric for
    the given EST index. </li>

    </ol>  
*/
typedef std::pair<int, float>  SMEntry;

/** \def SMList

    \brief Typedef for std::vector<SMEntry>.

    This typedef is a convenience definition to handle a vector of
    SMEntry objects associated with a given EST in the MST Cache.
*/
typedef std::vector<SMEntry>   SMList;

/** \def MSTCacheEntry

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

/** A cache to maintain similarity/distance metrics to build a Minimum
    Spanning Tree (MST).

    <p>This class provides some convenience methods to manage a cache
    of similarity metrics to facilitate rapid construction of a
    Minimum Spanning Tree (MST).  Furthermore, this class streamlines
    the operation of the Master and Worker processes when building a
    MST in a distributed manner.  The MSTCache caches similarity
    metrics to minimize the number of times the expensive operation of
    computing similarity metrics is performed via an EST
    analyzer. </p>

    \tparam<Comparator> The Comparator is a functor that is compatible
    with a binary predicate (or an adaptable binary predicate) that
    can be used to compare two floating point number to choose one of
    the them.  Examples are: std::less_equal<float> and
    std::greater_equal<float>.  In general if \c Comparator(x, y)
    returns \c true, then \c x is the result otherwise \c y is the
    result of comparison.  This comparator provides the MSTCache the
    feature to use either distance metrics (where smaller values are
    better) or similarity metrics (where larger values are better)
    without impacting the overall organization and functionality of
    this class.
*/
class MSTCache {
public:
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
    */
    MSTCache(const int totalESTCount, const int startOwnESTidx,
             const int numOwnESTs, const ESTAnalyzer *analyzer,
             const bool repopulateCaches);

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

    /** Clear caches that contain entries for a newly added estIdx.

        This method must be used to prune entries corresponding to the
        estIdxAdded from the various lists maintained by this cache. 
        The pruning operation is necessary for the following two
        purposes:

        <ol>

        <li>First it removes entries from various lists corresponding
        to estIdxAdded (because these entries are vestigial as the
        node as already been added to the MST) thereby reducing memory
        usage. </li>

        <li>Second, if a list become empty then the cache needs to be
        repopulated with fresh entries for further MST computations. </li>

        </ol>

        \param[in] estIdxAdded The index of the EST that has just been
        added to the MST.

        \param[out] repopulateList The index values of ESTs whose
        cache needs to be repopulated.  The list of EST indexes added
        to this list are only the ESTs that are owned by this process.
        If none of the caches need to be recomputed then this list is
        empty (all existing entries are removed).

	\param[in] prefixListSize If this parameter is set to true,
	then the number of entries added to the repopulateList is
	inserted at the beginning of the vector.  This eases
	transmission of lists via MPI to the master process.
    */
    void pruneCaches(const int estIdxAdded, std::vector<int>& repopulateList,
		     const bool prefixListSize = true);

    /** Merges a given SMList with current entries in the cache for a
        given EST.

        This method merges a given list with existing cache entries
        for the given EST such that only MaxCacheSize entries are
        retained in the cache.  The entries are chosen so that they
        have maximum similarity metrics.  Entries in the cache are
        removed to accommodate new entries.  Care is taken to ensure
        that the cache continues to remain sorted with the most
        similar entries at the top of the cache.

        \param[in] estIdx The index of the EST with which the given
        list must be merged.


        \param[in] list The SMList to be merged with any existing
        entries in the cache for the given estIdx.

        \param[in] maxCacheSize The maximum number of entries that
        must be retained in the cache after the list is merged.
    */
    void mergeList(const int estIdx, const SMList& list,
                   const int maxCacheSize);

    /** Obtains the top-most similar entry from the MSTCache.

        This method searches the similarity metrics stored in this
        cache for various ESTs to locate the best entry with the
        highest similarity.  It populates the parameters with the
        appropriate value.  Note that the parameters are initialized
        to -1, -1, and -1.0f respectively.

        \param[out] srcESTidx The source EST index from where the
        similarity metric is being measured.  The srcESTidx is already
        present in the MST.

        \param[out] destESTidx The destination EST index that is the
        best choice to be added to the MST (based on the local
        information).

        \param[out] metric The similarity/distance metric between the
        srcESTidx and the destESTidx.

    */
    void getBestEntry(int& srcESTidx, int& destESTidx, float& metric) const;

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
        similarity metric associated with each SMEntry in the list.

	\param[in] cacheSize The maximum size of this list after it
	has been sorted.
    */
    void sortAndPrune(SMList& list, const int cacheSize);

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

    /** Functor for SMEntry sorting.

        This Functor is used when sorting ESTs based on similarity
        metric by the sort method defined in this method.
    */
    class LessSMEntry : public std::binary_function<SMEntry, SMEntry, bool> {
    public:
        /** Constructor.

            The constructor requires a pointer to the ESTAnalyzer that
            is being used for analysis.  The analyzer is used to
            compare the metrics for sorting SMEntry objects.  This
            constructor is called in the MSTCache class for sorting
            cache entries.

            \param[in] analyzer The analyzer to be used for comparing
            the metric values associated with two SMEntry objects.
        */
        LessSMEntry(const ESTAnalyzer *analyzer) :
            comparator(analyzer) {}
        /** \fn operator()
            
            \brief operator<() for SMEntry.
            
            The following operator provides a convenient mechanism for
            comparing SMEntry objects for sorting.  This operator
            overrides the default STL operator() by comparing only the
            similarity metrics of two SMEntry objects (ignoring the
            EST indexes).  This method uses the metric comparison
            functionality provided by all EST analyzers.
            
            \param[in] sme1 The first SMEntry to be used for comparison.
            
            \param[in] sme2 The second SMEntry to be used for
            comparison.
            
            \return This method returns \c true if sme1.second is
            better than sme2.second.  Otherwise it returns \c false.
        */
        inline bool operator()(const SMEntry& sme1, const SMEntry& sme2)
        { return comparator->compareMetrics(sme1.second, sme2.second); }
        
    private:
        /** The functor for comparing.

            This comparator is set based on the template parameter
            associated with the enclosing class.
        */
        const ESTAnalyzer *const comparator;
    };

    static void copy_n(const SMList& input, const size_t count, 
                       SMList &output);

private:
    /**  The analyzer to be used for any EST comparison.
         
         This pointer is initialized in the constructor to point to
         the EST analyzer that is used for EST comparisons.  This
         pointer is essentially used for comparing EST metrics for
         ordering information in the cache.
    */
    const ESTAnalyzer *analyzer;
    
    /** The actual vector of cache entries for various nodes.

        This cache contains the list of cache entries for each node
        owned and managed by this MSTCache.  The list is initialized
        to have a set of empty lists in the constructor, when the
        cache is instantiated for use.
    */
    std::vector<MSTCacheEntry> cache;
    
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
        that is using this cache.  This value is used to internally
        normalize est index values to 0 for the first EST owned by
        this process.  This value is initialized in the constructor
        and is never changed during the lifetime of this object.
    */
    const int startOwnESTidx;

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
};

#endif
