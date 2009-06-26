#ifndef TRANS_MST_CLUSTER_MAKER_H
#define TRANS_MST_CLUSTER_MAKER_H

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

#include "MSTClusterMaker.h"
#include "TransCacheEntry.h"

/** \def HashMap<int, TransCacheEntry> MetricCacheMap

    \brief A shortcut to refer to a hash map of TransCacheEntry
    objects.

    This typedef provides a convenient shortcut to refer to a hash map
    containing the information to compute metrics using
    conditional-transitivity.  The key into this hash map is the index
    of the EST to which a CachedESTInfo corresponds to.  This index is
    the value in CachedESTInfo::estIdx.
*/
typedef HashMap<int, TransCacheEntry> MetricCacheMap;

/** A Minimum Spanning Tree (MST) based parallel cluster maker that
    uses conditional-transitivity relations to accelerate clustering.

    This class aims to enhance the performance (without \i
    signficiantly impacting quality of clusters) of the standard
    MSTClusterMaker by minimizing the number of the heavy weight
    analysis (such as: d2 or clu) performed for building the MST.
    Reduction in the number of required heavy weight analyses is
    achieved by applying the concept of transitivity in a conditional
    manner as follows:

    <ol>

    <li>Given three ESTs, say e<sub>1</sub>, e<sub>2</sub>, and
    e<sub>2</sub>, assume heavy weight metrics m<sub>12</sub> (between
    e<sub>1</sub> &amp; e<sub>2</sub> and m<sub>23</sub> (between
    e<sub>2</sub> &amp; e<sub>3</sub>) have been computed.</li>

    <li>When the heavy weight metric m<sub>13</sub> (between
    e<sub>1</sub> &amp; e<sub>3</sub>) is required, this class obtains
    this information in the following manner:

    <ol>

    <li>First a heuristic (such as: the <i>t/v</i> heuiristic) is
    applied to determine if the pair e<sub>1</sub> &amp; e<sub>3</sub>
    meet the stipulated relationship condition.  If not the ESTs are
    unrelated and conditional-transitivity is not
    applicable. Therefore, an invalid large distance value is assumed.</li>

    <li>However, if the pair of ESTs pass the heuristic test then they
    are related and conditional-transitivity is applicable.  In this
    case the heavy weight metric is deduced using a binary function
    <i>f</i>(m<sub>12</sub>, m<sub>23</sub>). Currently, the default
    implementation for <i>f</i> is \c std::max() method that selects
    the maximum of m<sub>12</sub> and m<sub>23</sub></li>
    
    </ol>
    
    </li>

    </ol>

    In order to rapidly retrieve pertinent metrics (or identify the
    lack of existing metrics), this class maintains a distributed set
    of hash maps in the \c transCache instance variable.  Each entry
    in the \c transCache hash map is a TransCacheEntry object and
    corresponds to a given EST (identified by its index).  The
    TransCacheEntry in-turn has a hash map of CachedESTInfo objects
    that provide all the necessary information to determine metrics
    via transitivity.
*/
class TransMSTClusterMaker : public MSTClusterMaker {
    friend class ClusterMakerFactory;
public:
    /** Method to display performance statistics.
        
        This method overrides the empty implementation in the base
        class to display additional statistics on
        conditional-transitivity application.

        \note This method calls the base class implementation
        first. Consequently, all the statistics displayed by the base
        class will continue to be displayed.
        
        \param[out] os The output stream to which the statistics must
        be written.
    */
    virtual void displayStats(std::ostream& os);
    
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~TransMSTClusterMaker();

protected:
    /** Helper method to call the actual heavy-weight analysis
        method(s).

        This method overrides the default implementation in the
        base-class before the actual call to the heavy weight analyzer
        is made.  This method attempts to use conditional-transitivity
        to try and minimize the number of calls to the heavy weight
        analyzer.  This method operates in the following manner:

        \param[in] otherEST The index of the other EST to which the
        metric is required.

        \return This method returns a similarity/distance metric by
        comparing the ESTs. This method may return -1, if the otherEST
        is significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
    */
    virtual float analyze(const int otherEST);
    
    /** Computes sends/receives similarity list for a given EST.

        This method overrides the default implementation in the base
        class that performs the core tasks of: (i) computing
        similarity metrics, (ii) gathering metrics from various
        processes to obtain the highest set of similarity metrics, and
        (iii) indicating to the manager that similarity metric
        calculations have been completed.

        \note This method is invoked on all the worker processes and
        the manager process.

        This method overrides the default implementation and performs
        the following tasks:

        <ol>

        <li>It calls the base class method to perform all the tasks in
        the usual manner.</li>

        <li>Next the process that owns \c estIdx broadcasts the cached
        metrics to all the other processes.  Each process in-turn
        waits to receive the cached metrics from the owner
        process.</li>

        <li>Then each process extracts and populates metrics for
        conditional-transitivity calculation by adding/updating
        appropriate entries to the transCache.
        
        </ol>
        
        \param[in] estIdx The index of the EST that was just added to
        the MST and for which the adjacent neighbors need to be
        determined.
        
        \param[out] smList If this pointer is not NULL, then this
        vector is populated with the set of metrics that were computed
        for estIdx <b>only on the owner process</b>.  This list
        contains the metrics collated from all the processes
        participating in the distributed computing process. Currently,
        this feature is used by TransMSTClusterMaker to obtain the
        list of metrics computed.        
    */
    virtual void populateCache(const int estIdx, SMList* metricList = NULL);
    
    /** The default constructor.
        
        The default constructor for this class.  The constructor is
        made protected so that this class cannot be directly
        instantiated.  However, since the ClusterMakerFactory is a
        friend of this class, an object can be instantiated via the
        ClusterMakerFactory::create() method.

        \param[inout] analyzer The EST analyzer to be used for
        obtaining similarity metrics between two ESTs.  This parameter
        is simply passed onto the base class.
        
	\param[in] refESTidx The reference EST index value to be used
	to root the spanning tree created by this method.  This
	parameter should be >= 0.  This value is simply passed onto
	the base class.

	\param[in] outputFile The name of the output file to which the
	raw MST cluster information is to be written.  If this
	parameter is the empty string then output is written to
	standard output.  This value is simply passed onto the base
	class.
    */
    TransMSTClusterMaker(ESTAnalyzer *analyzer, const int refESTidx,
                         const std::string& outputFile);

    /** Helper method to set the tvHeuristic pointer.

        This method is invoked from the populateCache() method to set
        a pointer to the tvHeuristic object in this class.  If the
        tvHeuristic pointer is not NULL, then this method returns
        immediately. In other words, the tvHeuristic pointer is set
        only once in this class, the first time this method is
        invoked.  If the tvHeuristic pointer is NULL, then this method
        operates as follows:

        <ol>

        <li>First, it searches the heuristicChain associated with the
        analyzer to determine if a "tv" heuristic class is present in
        the chain.  If so, the tvHeuristic pointer is made to point to
        that object.</li>

        <li>If the heuristicChain does not have a tv heuristic class,
        then this method creates a new object for further use. This
        object will eventually be deleted by the destructor.</li>

        </ol>
    */
    void setTVHeuristic();

    /** Helper method to remove invalid metrics from a given list of
	values.

	This is a helper method that is invoked from populateCache()
	method to remove invalid/unwanted entries from a given list of
	metrics. This method iterates over the entries in \c list and
	removes all entires that have their metric set to a value that
	is worse than badMetric value.

	\param[in] list The list of metrics that need to be
	pruned.

	\param[out] entries The list of entries that are better than
	badMetric.

	\param[in] badMetric The metric that defines the threshold and
	all entries that have an equal or worse value than badMetric
	are removed from the list.
    */
    void pruneMetricEntries(const SMList& list, SMList& entries,
			    const float badMetric);

    void processMetricList(SMList& remoteList);
    
private:
    /** An hash map that acts as the cache to hold known metrics.

        <p>This instance variable holds a cache of the previously
        computed metrics in a from that enables rapid application of
        conditional-transitivity to deduce metrics for related ESTs.
        The hash key for this hash map is the index of the reference
        EST for which metrics are to be computed.  That is, given two
        ESTs e<sub>i</sub> and e<sub>j</sub>,
        metricCache[i].getMetric(j, metric) provides an estimate of
        metric assuming e<sub>i</sub> and e<sub>j</sub> pass
        <i>t/v</i> heuristic.<p>

        <p>Entries in this cache are primarily added by the
        populateCache() method and entries are used by the analyze()
        method.</p>
    */
    MetricCacheMap metricCache;

    /** Instance variable to track total number of calls for analysis.

        This instance variable tracks the number of times the
        analyze() method was invoked. This instance variable
        essentially tracks the number of times the heavy-weight
        analysis method would have been invoked if
        conditional-transitivity was not applied. This instance
        variable is set to zero in the constructor, incremented in the
        analyze() method, and displayed by the displayStats() method.
    */
    int analyzeCount;
    
    /** Instance variable to track number of times transitivity was
        successfully applied.

        This instance variable tracks the number of times transitivity
        was successfully applied.  This number indices how effective
        transitivity approach was for a given set of ESTs. This
        instance variable is set to zero in the constructor,
        incremented in the analyze() method, and displayed by the
        displayStats() method.
    */
    int successCount;

    /** A pointer to the <i>t/v</i> heuristic object.

        This pointer holds a reference to a <i>t/v</i> heuristic
        object that used to establish conditional-transitivity between
        two given pairs of ESTs.  This object is initialized to NULL
        in the constructor and set to point to a valid tvHeuristic
        class in the setTVHeuristic() method (that is invoked from the
        populateCahce() method).  This pointer either points to a
        valid object in the heuristic chain associated with the
        analyzer or a custom heuristic object is created for use by
        this class.
    */
    Heuristic *tvHeuristic;

    /** The current reference EST for which analysis is underway.

        This instance variable tracks the index of the reference EST
        for which we are trying to figure out the closest neighbor in
        the MST.  This instance variable is set in the populateCache()
        method and is used in the analyze method. Note that the
        analyze method is actually called from within
        MSTClusterMaker::populateCache() (the base class
        implementation).
    */
    int currRefESTidx;
};

#endif
