#ifndef TRANS_MST_CLUSTER_MAKER_H
#define TRANS_MST_CLUSTER_MAKER_H

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
#include "TransCacheEntry.h"

#define NO_ERROR           0
#define ERROR_NO_HEURISTIC 1

/** \typedef std::vector<int, TransCacheEntry> MetricCacheMap

    \brief A shortcut to refer to a vector of TransCacheEntry objects.

    This typedef provides a convenient shortcut to refer to a vector
    containing the information to compute metrics using
    conditional-transitivity.  The index for this vector is the index
    of the EST to which a CachedESTInfo corresponds to.  This index is
    the value in CachedESTInfo::estIdx.
*/
typedef std::vector<TransCacheEntry*> MetricCacheMap;

/** A Minimum Spanning Tree (MST) based parallel cluster maker that
    uses conditional-transitivity relations to accelerate clustering.

    This class aims to enhance the performance (without \em
    signficiantly impacting quality of clusters) of the standard
    MSTClusterMaker by minimizing the number of the heavy weight
    analysis (such as: d2 or clu) performed for building the MST.
    Reduction in the number of required heavy weight analyses is
    achieved by applying the concept of transitivity in a conditional
    manner as follows:

    <ol>

    <li>Given three ESTs, say e<sub>1</sub>, e<sub>2</sub>, and
    e<sub>3</sub>, assume heavy weight metrics m<sub>12</sub> (between
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
    /** A method to handle initialization tasks for the TransMSTClusterMaker.

        This method is called after the ESTs have been loaded into the
        ESTAnalyzer and is used to initialize the vector of TransCacheEntrys,
        as well as to check for the existence of a heuristic chain
        on the analyzer.  Without a suitable heuristic, the
        TransMSTClusterMaker cannot function properly and will exit
        with an error code.
    */
    virtual bool initialize();
    
    /** Helper method to call the actual heavy-weight analysis
        method(s).

        This method overrides the default implementation in the
        base-class before the actual call to the heavy weight analyzer
        is made.  This method attempts to use conditional-transitivity
        to try and minimize the number of calls to the heavy weight
        analyzer.  This method operates in the following manner:

        \param[in] otherEST Pointer to the other EST to which the
        metric is required.

        \return This method returns a similarity/distance metric by
        comparing the ESTs. This method may return -1, if the otherEST
        is significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
    */
    virtual float analyze(const EST* otherEST);
    
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
        
        \param[out] metricList If this pointer is not NULL, then this
        vector is populated with the set of metrics that were computed
        for estIdx <b>only on the owner process</b>.  This list
        contains the metrics collated from all the processes
        participating in the distributed computing process. Currently,
        this feature is used by TransMSTClusterMaker to obtain the
        list of metrics computed.
    */
    virtual void populateCache(const int estIdx,
                               SMList* UNREFERENCED_PARAMETER(metricList)
                               = NULL);
    
    /** The default constructor.
        
        The default constructor for this class.  The constructor is
        made protected so that this class cannot be directly
        instantiated.  However, since the ClusterMakerFactory is a
        friend of this class, an object can be instantiated via the
        ClusterMakerFactory::create() method.

        \param[in,out] analyzer The EST analyzer to be used for
        obtaining similarity metrics between two ESTs.  This parameter
        is simply passed onto the base class.
    */
    TransMSTClusterMaker(ESTAnalyzer *analyzer);

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
    */
    void pruneMetricEntries(const SMList& list, SMList& entries);

    /** Helper method to process a list of similarity metrics for caching.
        This method is called in the populateCache method and will take
        care of a list of similarity metrics, either computed locally or
        received via remote communication.  This method goes through the list
        of similarity metrics and updates TransCacheEntrys as appropriate.
        
        \param[in] metricList The list of metrics received or computed.
    */
    void processMetricList(SMList& metricList);
    
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

    float badMetric;
};

#endif
