#ifndef MST_CLUSTER_MAKER_H
#define MST_CLUSTER_MAKER_H

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

#include "ClusterMaker.h"
#include "MSTCache.h"
#include "MSTCluster.h"
#include "MST.h"

#include <fstream>

/** A constant to remove magic 0 (zero) success checks. */
const int NO_ERROR = 0;

/** Base class for all Minimum Spanning Tree (MST) based parallel
    cluster makers.

    This class serves as a generic base class for the different
    cluster makers to construct a MST-based EST clusters in a
    parallel/distributed manner using the Message Passing Interface
    (MPI) library.  This class includes functionality for both the
    Manager (MPI Rank == 0) and Worker (MPI Rank > 0) processes.
    Necessary functionality to distinguish and operate either as
    Manager or Worker is already built into the class.  This class
    uses the MSTCache and MSTCluster classes to help in performing the
    various activities.  Refer to the documentation on the various
    method for detailed description on their functionality and
    usage.
*/
class MSTClusterMaker : public ClusterMaker {
public:
    /** The set of tags exchanged between various processes.

        This enum provides meanigful names to the various tags
        (integers) exchanged between the master and worker processes
        participating in the construction of a MST in a
        parallel/distributed manner.
    */
    enum MessageTags{REPOPULATE_REQUEST, COMPUTE_SIMILARITY_REQUEST,
                     SIMILARITY_LIST, SIMILARITY_COMPUTATION_DONE,
                     COMPUTE_MAX_SIMILARITY_REQUEST, MAX_SIMILARITY_RESPONSE,
                     ADD_EST, TRANSITIVITY_LIST, COMPUTE_TOTAL_ANALYSIS_COUNT,
                     CLUSTER_INFO_LIST};
    
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~MSTClusterMaker();

   /** Add valid command line arguments for this analyzer.

        This method must be used to add all valid command line options
        that are supported by this analyzer.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.

        \note Derived EST analyzer classes may override this method to
        display help for their custom command line arguments.  When
        this method is overridden don't forget to call the
        corresponding base class implementation to add common options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

    /** A method to handle initialization tasks for the MSTClusterMaker.
		
        This method is called after the ESTs have been loaded and the
        clustering is just about to commence. This method calls the
        corresponding method in the base class to ensure base class
        initialization proceeds successfully (call to base class
        inititalizes the ESTAnalyzer and sets up estList pointer).  It
        then validates command-line parameters.

		\return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    virtual bool initialize();
	
    /** Method to begin clustering.

        This method must be used to create clusters based on a given
        EST analysis method.  This method performs the following
        tasks:
        
    */
    virtual int makeClusters();

    /** Method to display performance statistics.

        This method overrides the empty implementation in the base
        class to display statistics on cache usage and MPI calls for
        tracking and reporting the performance and behavior of this
        class.

        \param[out] os The output stream to which the statistics must
        be written.        
    */
    virtual void displayStats(std::ostream& os);

    /** Add a dummy cluster to the cluster maker.
        
        This method can be used to add a dummy cluster to the cluster
        maker. The dummy clusters are added as direct descendants of
        the \c root cluster with the given name.

        \note This method is currently used by the Filter hierarchy to
        add ESTs that are logically filtered out and must not be part
        of the core clustering process.
        
        \param[in] name A human readable name to be set for this
        cluster. No special checks are made on the contents of the
        string.
        
        \return This method returns the ID value set for the newly
        added dummy index.
    */
    virtual int addDummyCluster(const std::string name);

    /** Add a EST directly to a given cluster.

        This method can be used to add an EST directly to a
        cluster. This bypasses any traditional mechanism and directly
        adds the EST to the specified cluster.

        \note The EST is added with an invalid metric value. ESTs
        added to a cluster are not included in the standard clustering
        process. Adding an EST that has already been added to the
        same/another cluster results in undefined behaviors.

        \param[in] clusterID The unique ID of the cluster to which the
        EST is to be added. This value must have been obtained from an
        earlier (successful) call to the ClusterMaker::addDummyCluster
        method.
        
        \param[in] estIdx The EST to be added to the given
        cluster. Once the EST has been added to this cluster it will
        not be included in the clustering process performed by this
        cluster maker.
    */
    virtual void addEST(const int clusterID, const int estIdx);

    /** Obtain a pointer to the root cluster.

        This method can be used to obtain a reference to the set of
        clusters that has been built by this class.  This method
        returns a valid set of clusters only after the makeClusters()
        method has successfully completed.

        \return A pointer to the root of the clusters built by this
        cluster maker.
    */
    const MSTCluster* getClusters() const { return &root; }

    /** Obtain a mutable-reference to the root cluster.

        This method can be used to obtain a reference to the set of
        clusters that has been built by this class.  This method
        returns a valid set of clusters only after the makeClusters()
        method has successfully completed.

        \return A reference to the root of the clusters built by this
        cluster maker.
    */
    MSTCluster& getClusters() { return root; }
    
    /** Obtain a pointer to the Minimum Spanning Tree (MST) built by
        this class.

        This method implements the virtual method in the base class.
        It can be used to obtain an MST object built by this class
        (with help from one of the derived child classes).

        \return The MST object built by this class.  The returned
        pointer can be NULL if an MST has not been built.
    */
    const MST* getMST() const { return mst; }
	
protected:
    /** Variable to indicate per-EST similarity cache size.

        This variable is used to indicate the number of similarity
        metrics that must be cached for a given EST.  This value is
        initialized to 128.  The value can be changed by the user via
        the \c --cache command-line argument.
    */
    int cacheSize;

    /** Variable to indicate if strict ordering of worker Ranks must
        be followed.

        <p> If this member variable is \c true, then messages
        dispatched by workers and the manager are always read in a
        fixed order of increasing ranks.  That is, messages from rank
        0 (zero) are processed first, then messages from process with
        rank 1, so on and so forth.  On the other hand if this
        variable is \c false, then messages are processed in the order
        they are received. </p>
        
        <p> The strictOrder approach guarantees consistent results for
        each run (involving the same number of processes) and the
        resulting MSTs are alll identical.  However, a process may
        have to wait (idle wasting time) until a message from the
        appropriate process (with a given rank) is actually received.
        This may slow down the overall computational rate,
        particularly when the work load get's skewed toward the end of
        MST construction.</p>

        <p> On the other hand, if strictOrder is relaxed (by setting
        strictOrder variable to \c false) then messsages are processed
        as soon as they are received, in the order in which messages
        arrive.  This approach minimizes wait times.  However, the MST
        constructed between multiple runs may not be identical as
        equidistant (or nodes with same similarity metrics) nodes may
        be processed in different order.  Reordering of equidistant
        nodes occur because in this mode a total order is not enfored
        and only a partial order of nodes is performed.</p>

        <p> By default strictOrder is enabled.  However, the value can
        be changed by the user through the command line argument \c
        --no-order.</p>
    */
    bool strictOrder;

    /** Command line option to avoid the clustering phase.

        If this member variable is \c true, then this class only
        generates MST information and does not do clustering.  By
        default this variable is initialized to \c false.  However,
        the value can be changed by the user through the command line
        argument \c --dont-cluster.
    */
    bool dontCluster;

    /** Variable to indicate if MST information must be simply read
        from a given file.

        This member variable is used to hold the name of the file
        (with full path) from where MST information must be read.
        This instance variable is initialized to an empty string.
        However, if the input MST file is specified (via the
        command-line argument \c --input-mst-file) then MST building
        is skipped and MST data read from the file is used for further
        processing.
    */
    std::string inputMSTFile;

    /** Command line option to suppress cache repopulation.
	
        <p> If this member variable is \c true, then this class does
        not repopulate caches once a EST cache becomes empty.  By
        default this variable is initialized to \c false.  However,
        the value can be changed by the user through command line
        argument \c --no-cache-repop.</p>

        <p> If this parameter is not specified then the MSTCache will
        request lists to be repopulated when needed.  Repopulating
        lists guarantees that ultimately a MST will be developed.  If
        repopulation is suppressed via this parameter then the
        resulting spanning tree may not be a MST; however computation
        time decreases. </p>
    */
    bool noCacheRepop;

    /** Command line option to enable maximum use of precomputed
        scores for building MST.

        <p>If this member variable is set to a value other than -1,
        then the MSTClusterMaker will try to use all the ESTs that
        have a metric better than a given threshold specified in
        maxUse. Maximally using good metrics will ultimately reduce
        the total number of analysis that need to be performed,
        thereby reducing overall time for clustering.  This value can
        be changed via the command-line argument \c --maxUse \e
        threshold.
    */
    float maxUse;
    
    /** Command line option to set the clustering threshold to be
        used in deriving clusters from the MST.
        
        This member variable is used to set the clustering threshold.
        This value can be set via the command-line parameter \c
        --clsThreshold.  There are two "special" options here:
        
        <ul>

        <li>\b 1.0 : Corresponds to the TwoPassD2 analyzer which uses
        different window lengths, each with different thresholds.  If
        TwoPassD2 is used then \c --clsThreshold clsThreshold must be
        set to 1.0.</p>
        
        <li> \b -1 : Corresponds to the mean/variance-based threshold.
        Intended to be used with the CLU analyzer.</li>

        </ul>
    */
    float clsThreshold;

    /** Command line option to set the type of cache to be used by
        PEACE.

        This member variable is used to indicate the type of cache
        that must be used to store metrics to facilitate rapid
        construction of the MST.  The default cache used in the
        MSTHashCache indicated by the cacheType set to \c "hash".  The
        alternative cache in the MSTMultiListCache (indicated by
        cacheType value of \c "mlist").  The user may override the
        default using the command line parameter \c --cacheType.
    */
    std::string cacheType;

	/** Name of file to report progress in during MST construction.

        This command line argument provides the name of the log file
		where progress information is to be written. The progress
		information is in the form: \#estsProcessed, \#ests. This
		value is specified via the command line argument \c
		--progress.  This feature is used by PEACE GUI.
    */
    std::string progFileName;
	
    /** Helper method to perform manager tasks.

        This method has been introduced to streamline the operations
        of the MSTClusterMaker when it operates as the manager. The
        MPI process with Rank 0 (zero) acts as the manager and
        coordinates all the activities of the MSTClusterMaker.  This
        method is invoked from the makeClusters() method.

        \return This method returns 0 (zero) if clusters were created
        successfully.  Otherwise this method returns a non-zero value
        indicating an error.
    */
    virtual int manager();

    /** Helper method to perform worker tasks.

        This method has been introduced to streamline the operations
        of the MSTClusterMaker when it operates as a worker. All the
        MPI processes with non-zero rank act as a worker and
        collaborate with the manager to assist in various activities
        of the MSTClusterMaker.  This method is invoked from the
        makeClusters() method.

        \return This method returns 0 (zero) if clusters were created
        successfully.  Otherwise this method returns a non-zero value
        indicating an error.
    */
    virtual int worker();

    /** Helper method to generate (or compute) or load MST data from file.

        This is a helper method that is invoked from the
        makeClusters() method. This method performs one of the
        following tasks:

        <ul>

        <li>If an #inputMSTFile has not been specified, then this
        method builds an MST (either on one process or many MPI
        processes).  This method first creates a local cache local
        cache that contains information to build the MST. It then
        builds the MST calling the manager() or worker() method
        depending on the MPI-rank of this process.</li>

        <li>If an #inputMSTFile has indeed been specified as a command
        line parameter, then this method loads the MST from the
        specified MST file.</li>
        
        </ul>

        \return This method returns zero on success indicating that a
        valid MST is available for further processing. Otherwise this
        method returns a non-zero code signalling error.
    */
    virtual int populateMST();

    /** Utility method to do the final clustering step.

        This is a refactored (primarily to keep the code clutter to a
        minimum) utility method that is used to perform the final step
        in clustering.  This method essentially calls the
        MSTCluster::makeClusters method that builds the clusters using
        the MST.
        
        \return This method returns zero on success. If errors occur,
        this method returns a non-zero error code.
    */
    virtual int buildClusters();
    
    /** Helper method to call the actual heavy-weight analysis
        method(s).

        This is a helper method that is invoked from the
        populateCache() method to obtain the relationship metric
        (either via CLU or d2) between the current parent EST and the
        given otherEST. This method was introduced to enable chlid
        classes (such as TransMSTClusterMaker) to conveniently
        intercept analyzer calls and potentially shortcircuit them
        using concepts of conditional-transitivity.

        \param[in] otherEST An immutable pointer to the other EST
        object to which the metric is required.

        \return This method returns a similarity/distance metric by
        comparing the ESTs. This method may return -1, if the otherEST
        is significantly different from the reference EST (possibly
        warranting no further analysis) that a meanigful metric cannot
        be generated.
    */
    virtual float analyze(const EST *otherEST);

    /** Helper method to compute \i local distance/similarity metrics
        and populate the information in the supplied smList.

        This method is invoked from the populateCache() method in this
        class. This method has been introduced to streamline the code
        in populateCache() method and provide a mechanism for derived
        classes to further customize its operations.  The
        implementation for this method in this class operates as
        follows:

        <ol>

        <li>Each process computes the range of local ESTs to operate
        on via a call to the getLocallyOwnedESTidx() method.</li>

        <li>It set's up the supplied estIdx as the reference EST for
        the series of analyses being performed in the next step.</li>
        
        <li>For each locally-owned \i unprocessed EST that is not in
        the MST, this method:

        <ol>

        <li> Computes the similarity/distance metric via the analyzer.</li>

        <li>Furthermore, any alignment information and direction data
        supplied by the analyzer (if any) is obtained and stored into
        a CachedESTInfo object.</li>

        <li>The computed CachedESTInfo object is then added to the
        supplied smList (second parameter). The method ensures that
        the list has at least one entry in it if an analysis was
        performed.</li>

        </ol>
        
        </ol>

        \param[in] estIdx The index of the reference EST (typically
        the EST that was just added to the MST) for which the adjacent
        neighbors (based on similarity/distance metric) need to be
        determined.

        \param[out] smList The vector to which neighbors are to be
        stored and returned.  This list is not cleared by this
        method. Moreover, the list is not sorted or organized in any
        manner.
    */
    virtual void computeSMList(const int estIdx, SMList& smList);
    
    /** Computes sends/receives similarity list for a given EST.

        This method is a shared method that is used by both the
        manager and workers.  This method is used to compute the
        similarity metric and cache the highest set of similarity
        metrics.  This method operates as follows:

        <ol>

        <li>Each process computes a subset of the EST similarity
        metric in the range \em k*Rank < otherEstIdx < (\em k+1)*Rank,
        where k=estList.size() / MPI::COMM_WORLD.Get_size(), and Rank
        is the MPI rank of this process. </li>

        <li>If this process is the cache owner for the est, (that is,
        estIdx % Rank == 0), then it receives data from other
        processes and merges the information with its own list,
        retaining the top-most similarity metrics.</li>

        <li>If this process is \b not the cache owner for the est,
        (that is, estIdx % Rank != 0), then it sends the computed
        similarity metrics to the owner process. <li>
        
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
    virtual void populateCache(const int estIdx, SMList* metricList = NULL);

    /** Method to generate progress logs (if requested by user).

        This method is a helper method that is called from the core
        manager() method loop to generate progress logs as ESTs are
        analyzed and updated.  This method cuts logs only if the
        progFileName comamnd line argument was specified and the
        progress file could be created.

        \note The progress file is opened only once, first time this
        method is called. So this method should not hurt performance
        by much.
        
        \param[in] estsAnalyzed The number of ESTs analyzed thus far.

        \param[in] totalESTcount The total number of ESTs to be
        analyzed.

        \param[in] estAdded Index of the EST that was just added (if
        any).

        \param[in] metric The metric associated with the entry just
        added.

        \param[in] direction The direction flag (1: FC, -1: RC) for
        the entry just added.
    */
    void updateProgress(const int estsAnalyzed,   const int totalESTcount,
                        const int estAdded  = -1, const float metric = -1,
                        const int direction = 0);
        
    /** Helper method in Manager process to update distributed caches.

        This is a helper method that is used only in the Manager
        process to perform the following tasks using the newly added
        estIdx value:

        <ol>

        <li>First, this method broadcasts the newly added EST index
        (\c estIdx) to all the workers.</li>

        <li>Next it prunes it local cache via the
        MSTCache::pruneCaches() method.</li>

        <li>It then collects requests to repopulate specific caches
        from all the workers.</li>

        <li>It then adds the newly created est to the list of caches
        to be repopulated and broadcasts request to repopulate caches
        to each worker and participates in cache repopulation task by
        calling the populateCache() method.</li>

        </ol>

        \param[in] estIdx The index of the newly added EST.

        \param[in] refreshEST If this flag is \c true (the default
        value), then the neighbors for the newly added EST (specified
        by estIdx) are computed and the caches are updated.
        
        \return This method returns 0 on success or an suitable error
        code on failure.
    */
    int managerUpdateCaches(int estIdx, const bool refreshEST = true);

    /** Helper method in \b Manager process to collaboratively compute
		the next EST to be added to the MST.

        This is a helper method that is used only in the Manager
        process to perform the following tasks using the newly added
        estIdx value:

        <ol>

        <li>First, this method sends request to compute the best local
        choice to each of the worker processes.</li>

        <li>Next it computes its own local (at the Manager's end) best
        choice for the next EST node to be added.</li>

        <li>It then collects response for best local choice from each
        worker process and tracks the best reported value.</li>

        </ol>

        \param[out] parentESTidx The source EST index from where the
        similarity metric is being measured.  The srcESTidx is already
        present in the MST.

        \param[out] estToAdd The destination EST index that is the
        best choice to be added to the MST (based on the local
        information).

        \param[out] similarity The similarity metric between the
        srcESTidx and the destESTidx.
        
        \param[out] alignmentData The alignment information between
        the two ESTs represented by their index values in parentESTidx
        and estToAdd.

        \param[out] directionData The direction information between
        the two ESTs represented by their index values in parentESTidx
        and estToAdd.
    */
    void computeNextESTidx(int& parentESTidx, int& estToAdd,
                           float &similarity, int& alignmentData,
                           int& directionData) const;
    
    /** Determine the owner process Rank for a given estIdx.

        This method is a convenience method to determine the Rank of
        the process that logically owns a given EST.  The owning
        process is responsible for maintaining the cache for a given
        EST.  The owners are assigned in a simple fashion and ESTs are
        evenly divided up amongst all the processes.

        \param[in] estIdx The index of the EST whose owner process's
        rank is requested.  It is assumed that the estIdx is valid.
        If invalid EST index values are supplied then the operation of
        this method is undefined.

        \note This method must be invoked only after MPI::Intialize()
        has beeen called and the ESTs to be processed have be loaded
        (so that EST::getESTList() returns a valid list of ESTs).
        
        
        \return The rank of the owner process for the given estIdx.
    */
    int getOwnerProcess(const int estIdx) const;

    /** Helper method to compute the start and ending indexes of the
        EST that this process owns.

        This method was introduced to keep the math and logic clutter
        involved in computing the list of owned ESTs out of the
        methods that use the information.  This method returns the
        range, such that: \c startIndex <= \em ownedESTidx < \c
        endIndex.

        
        \note This method must be invoked only after MPI::Intialize()
        has beeen called and the ESTs to be processed have be loaded
        (so that EST::getESTList() returns a valid list of ESTs).

        \param[out] startIndex The starting (zero-based) index value
        of the contiguous range of ESTs that this process owns.

        \param[out] endIndex The ending (zero-based) index value of
        the contiguous range ESTs that this process owns.  The value
        returned in this parameter is \b not included in the range of
        values.
    */
    void getLocallyOwnedESTidx(int& startIndex, int& endIndex);

    /** Helper method for a worker process.

        This method is invoked from the worker() method to receive and
        process various requests from the manager process.  This
        method currently handles the following requests:
        
        <ul>
        
        <li>\c COMPUTE_SIMILARITY_REQUEST : Computes the subset of the
        similarity metric for the given EST index and returns the
        partial list back to the owner process.</li>

        <li> \c COMPUTE_MAX_SIMILARITY_REQUEST : Computes the highest
        similarity value between all the ESTs on this cluster and
        returns the top entry back to the manager.  Once this request
        has been processed this method returns control back.
        
        </ul>
    */
    void workerProcessRequests();

    /** Distribute data and tag to all the workers.

        This method provides a convenient mechanism to broadcast a
        given integer data and tag to all the workers.

        \param[in] data The integer to be sent to each and every
        worker.

        \param[in] tag The message tag to be sent to each and every
        process.
    */
    void sendToWorkers(int data, const int tag) const;

    /** Method to detect if a given SMList has at least one, valid
        entry.
        
        This method is used (in the populateCache()) to determine if a
        given SMList has at least one valid entry.  This method is
        useful particularly when a empty SMList is received from a
        remote process and in this case there ine one entry in the
        SMList (-1, -1). 
        
        \param[in] list The list to check if it has a valid entry.
    */
    bool hasValidSMEntry(const SMList& list) const;
    
    /** Helper method to distribute index of newly added EST to all
        workers and gather cache repopulation requests.

        This is a helper method that was added to streamline the code
        in managerUpdateCaches method.  This method performs the
        following tasks:

        <ol>

        <li>First it uses the \c sendToWorkers() method to distribute
        the \c estIdx (parameter) value to all the workers. </li>

        <li> Next it prunes the local caches on the manager.</li>

        <li>It then obtains repopulation requests from each worker and
        places EST indexes to be repopulated in the repoulateList
        parameter.</li>

        </ol>

        \note This method must be inovked only on the manager.

        \param[in] estIdx The index of the newly added EST that must
        be distributed to all the workers.

        \param[out] repopulateList A vector that will contain the list
        of ESTs that need to be repopulated (based on requests
        received from various workers).
    */
    void estAdded(const int estIdx, std::vector<int>& repopulateList);

    /** Helper method in \b Manager process to add as many child nodes
        as possible for the given parent.

        This is a helper method that is used only in the Manager
        process <b>only when the \c maxUse parameter is != -1</b>.
        This method tries to add more children rooted at the given
        parent to the MST as long as the metric is better than \c
        maxUse value.  This method operates as follows:

        <ol>

        <li>First, this method sends request to compute the best local
        choice to each of the worker processes.</li>

        <li>Next it computes its own local (at the Manager's end) best
        choice for the next EST node to be added.</li>

        <li>It then collects response for best local choice from each
        worker process and tracks the best reported value.</li>

        <li>If the next best entry is still rooted at this parent and
        the metric is better than \c maxUse then the EST is added to
        MST and the process is repeated from step 1. Otherwise, the
        parameters are updated to the last added EST and the method
        returns.</li>

        </ol>

        \param[in] parentESTidx The source EST index from where the
        similarity metric is being measured.  The parentESTidx is
        already present in the MST.

        \param[in,out] estToAdd The EST that has just been added to
        the MST.  This method updates this value if additional ESTs
        are added to the MST by this method.

        \param[in,out] metric The similarity/distance metric between
        the parentESTidx and the estToAdd.  This method updates this
        value if additional ESTs are added to the MST by this method.
        
        \param[in,out] alignmentData The alignment information between
        the two ESTs represented by their index values in parentESTidx
        and estToAdd.  This method updates this value if additional
        ESTs are added to the MST by this method.
        
        \param[in,out] directionData The direction information between
        the two ESTs represented by their index values in parentESTidx
        and estToAdd.  This method updates this value if additional
        ESTs are added to the MST by this method.
        
        \param[in,out] pendingESTs The number of pending ESTs that
        have not yet been added to the MST.  This value is used and
        udpated by this method each time it adds a EST.

        \param[in] TotalESTcount The total number of ESTs being
        processed.  This value is used to update progress/status
        information.
    */
    void addMoreChildESTs(const int parentESTidx, int& estToAdd,
                          float &metric, int& alignmentData,
                          int& directionData, int& pendingESTs,
                          const int TotalESTcount);
    
    /** The constructor.
        
        The only constructor for this class.  The constructor is made
        protected so that this class cannot be directly instantiated.
        Instead one of the derived classes must be instantiated
        and used.

        \param[in] name The name to be associated with this cluster
        maker.  This value is simply passed onto the base class and is
        used when generating messages.
        
        \param[in,out] analyzer The EST analyzer to be used for
        obtaining similarity metrics between two ESTs.  This parameter
        is simply passed onto the base class.
    */
    MSTClusterMaker(const std::string& name, ESTAnalyzer *analyzer);
    
    /** The cache that holds similarity metrics for MST construction.

        This object is used to cache the similarity metrics for all
        ESTs that are owned by this process (that is, estIdx % Rank ==
        0, where Rank is the MPI rank of this process).  The cache
        contains similarity metrics to facilitate rapid construction
        of the MST.  Both the manager and worker processes have their
        own caches and manage them independently.  This spreads out
        the memory requirement for the caches across multiple
        processes enabling large (in 10s of GB) caches.

        The cache is created just before the clustering process
        commences and is deleted immediately after the clustering
        process (to minimize memory footprint).
    */
    MSTCache *cache;

    /** File stream to log progress information.

        This output stream is created when the first progress information
        is logged and closed after the last progress information has been
        logged. The progress information is generated by the updateProgress
        method if progressFileName is not NULL.
    */
    std::ofstream progressFile;
    
    /** The output file to which the clustering results are to be
        written.

        This isntance variable tracks the output file to which the
        data is to be written.  This option should be moved to the
        output sub-system.
    */
    std::string outputFileName;

    /** Flag to indicate if progress information should also be
        printed on console.

        This flag is set via the \c --print-progress command-line
        argument.
    */
    bool printProgress;

private:    
    /** The Minimum Spanning Tree (MST) built by this class.

        This instance variable holds a pointer to the MST created by
        this class when it operates as a manager process.  This
        pointer is initialized to NULL and a MST is created in the
        manager() method.
    */
    MST* mst;

    /** The top-level root cluster that contains all other clusters.

        This member represents the top-level root cluster that contain
        all other clusters created by this cluster maker. This cluster
        also contains dummy clusters that are created by Filter
        objects used in conjunction with clustering.
    */
    MSTCluster root;

    /** The hint key that is used to add hint for normal or
        reverse-complement D2 computation.
	
        This hint key is used to obtain a hint for the \c MST_RC in
        the \c hints hash map. This string is defined as a constant to
        save compute time in the core \c runHeuristics method.
    */
    const std::string hintKey_MST_RC;
};

#endif
