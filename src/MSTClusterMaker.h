#ifndef MST_CLUSTER_MAKER_H
#define MST_CLUSTER_MAKER_H

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

#include "ClusterMaker.h"
#include "MSTCache.h"
#include "MST.h"

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
class MSTClusterMaker : public ClusterMaker {
    friend class ClusterMakerFactory;
public:
    /** The set of tags exchanged between various processes.

        This enum provides meanigful names to the various tags
        (integers) exchanged between the master and worker processes
        participating in the construction of a MST in a
        parallel/distributed manner.
    */
    enum MessageTags{REPOPULATE_REQUEST, COMPUTE_SIMILARITY_REQUEST,
		     SIMILARITY_LIST, SIMILARITY_COMPUTATION_DONE,
		     COMPUTE_MAX_SIMILARITY_REQUEST, MAX_SIMILARITY_RESPONSE};
    
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    ~MSTClusterMaker();

    /** Display valid command line arguments for this cluster maker.

        This method must be used to display all valid command line
        options that are supported by this cluster maker (and its base
        classes).
        
        \note This method calls the base class's showArguments first.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.

        This method is used to process command line arguments specific
        to this cluster maker.  This method is typically used from the
        main method just after the cluster maker has been
        instantiated.  This method consumes all valid command line
        arguments.  If the command line arguments were valid and
        successfully processed, then this method returns \c true.

        \note This method consumes its custom command line arguments
        first and then call's the base class's parseArguments() method.
        
        \param[inout] argc The number of command line arguments to be
        processed.

        \param[inout] argc The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.
    */
    virtual bool parseArguments(int& argc, char **argv);

    /** Method to begin clustering.

        This method must be used to create clusters based on a given
        EST analysis method.  This method performs the following
        tasks:
        
    */
    virtual int makeClusters();
    
protected:
    /** Variable to indicate per-EST similarity cache size.

        This variable is used to indicate the number of similarity
        metrics that must be cached for a given EST.  This value is
        initialized to 128.  The value is changed by the
        parseArguments() method if the user has specified an option to
        override the default.
    */
    static int cacheSize;

    /** Command line option to set percentile value to compute
	clustering threshold.

	<p> This variable is used to indicate the percentile value
        that must be used to determine the threshold for clustering.
        This value is initialized to 1.0. This value is ultimately
        used in the MSTCluster::calculateThreshold() method to compute
        the threshold using the formula: </p>

	<i> threshold = mean + (stDev * percentile); </li>

	<p> The value is changed by the parseArguments() method if the
        user has specified the --percentile option to override the
        default. </p>
    */
    static double percentile;
        
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
        be changed by the user through command line arguments.  The
        change of value occurs in the parseArguments() method if the
        user has specified an option to override the default.
    */
    static bool strictOrder;

    /** Command line option to avoid the clustering phase.

        If this member variable is \c true, then this class only
        generates MST information and does not do clustering.  By
        default this variable is initialized to \c false.  However,
        the value can be changed by the user through command line
        arguments.  The change of value occurs in the parseArguments()
        method if the user has specified an option to override the
        default.
    */
    static bool dontCluster;

    /** Command line option to print a pretty cluster tree.

        If this member variable is \c true, then this class prints a
        pretty ASCII tree with the cluster information By default this
        variable is initialized to \c false.  However, the value can
        be changed by the user through command line argument
        (--pretty-print).  The change of value occurs in the
        parseArguments() method if the user has specified an option to
        override the default.
    */
    static bool prettyPrint;
    
    /** Variable to indicate if MST information must be simply read
        from a given file.

        This member variable is used to hold the name of the file
        (with full path) from where MST information must be read.
        This instance variable is initialized to NULL.  However, if
        the input MST file is specified then MST building is skipped
        and MST data read from the file is used for further processing.
    */
    static char* inputMSTFile;

    /** Variable to indicate if MST information must be written to a
        given file.

        This member variable is used to hold the name of the file
        (with full path) to which MST information must be written.
        This instance variable is initialized to NULL.  However, if
        the output MST file is specified then MST data built by this
        program is written to the specified file.
    */
    static char* outputMSTFile;

    /** Command line option to suppress cache repopulation.
	
        <p> If this member variable is \c true, then this class does
        not repopulate caches once a EST cache becomes empty.  By
        default this variable is initialized to \c false.  However,
        the value can be changed by the user through command line
        argument (--no-cache-repop).  The change of value occurs in
        the parseArguments() method if the user has specified an
        option to override the default.</p>

	<p> If this parameter is not specified then the MSTCache will
	request lists to be repopulated when needed.  Repopulating
	lists guarantees that ultimately a MST will be developed.  If
	repopulation is suppressed via this parameter then the
	resulting spanning tree may not be a MST; however computation
	time decreases. </p>
    */
    static bool noCacheRepop;
    
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

    /** Computes sends/receives similarity list for a given EST.

        This method is a shared method that is used by both the
        manager and workers.  This method is used to compute the
        similarity metric and cache the highest set of similarity
        metrics.  This method operates as follows:

        <ol>

        <li>Each process computes a subset of the EST similarity
        metric in the range \i k*Rank < otherEstIdx < (\i k+1)*Rank,
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
    */
    virtual void populateCache(const int estIdx);

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

        <il>It then adds the newly created est to the list of caches
        to be repopulated and broadcasts request to repopulate caches
        to each worker and participates in cache repopulation task by
        calling the populateCache() method.</li>

        </ol>

        \param[in] estIdx The index of the newly added EST.

        \return This method returns MPI::SUCCESS on success or an
        suitable error code on failure.
    */
    int managerUpdateCaches(int estIdx);

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

        \param[out] srcESTidx The source EST index from where the
        similarity metric is being measured.  The srcESTidx is already
        present in the MST.

        \param[out] destESTidx The destination EST index that is the
        best choice to be added to the MST (based on the local
        information).

        \param[out] similarity The similarity metric between the
        srcESTidx and the destESTidx.
    */
    void computeNextESTidx(int& parentESTidx, int& estToAdd,
			   float &similarity) const;
    
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
        range, such that: \c startIndex <= \i ownedESTidx < \c
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
    void getOwnedESTidx(int& startIndex, int& endIndex);

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
    
private:
    /** The default constructor.

        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ClusterMakerFactory is a
        friend of this class, an object can be instantiated via teh
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
    MSTClusterMaker(ESTAnalyzer *analyzer, const int refESTidx,
                    const std::string& outputFile);

    /** The set of common arguments for the MST cluster maker.

        This instance variable contains a static list of arguments
        that are common all the MST cluster maker objects.
    */
    static arg_parser::arg_record argsList[];

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

    /** The Minimum Spanning Tree (MST) built by this class.

        This instance variable holds a pointer to the MST created by
        this class when it operates as a manager process.  This
        pointer is initialized to NULL and a MST is created in the
        manager() method.
    */
    MST* mst;
};

#endif