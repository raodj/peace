#ifndef NN_MST_CLUSTER_MAKER_H
#define NN_MST_CLUSTER_MAKER_H

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

/** A Nearest-Neighbor (NN) and Minimum Spanning Tree (MST) based
    parallel cluster maker.

    <p>This class encapsulates the core functionality needed to
    construct a MST-based EST clusters in a parallel/distributed
    manner using the Message Passing Interface (MPI) library.  The
    specific additional functionality provided by this class (that
    distinguishes it from the underlying default approach) is the the
    use of a nearest-neighbor pool of candidate ESTs.  Specifically,
    rather than inspecting all pertinent ESTs to determine best match,
    this class maintains a working pool of nearest-neighbor ESTs that
    serve as primary candidates to be checked to identify similar ESTs
    to be clustered togehter. The nearest-neighbor pool is managed and
    used in the following manner:

    <ul>

    <li>Whenever the nearest-neighbor pool is empty, all ESTs (that
    are not already in the MST) are implicitly assumed to be
    candidates and are processed.</li>

    <li>As candiate ESTs are processed, ESTs that have a score better
    than ESTAnalyzer::getInvalidMetric(), are added to the
    nearest-neighbor pool for consideration in the next round.</li>

    <li>If all entries in the explicit nearest-neighbor pool are
    analyzed without encountering at least one successful candiate,
    then remaining ESTs are considered as candiates, they are
    processed, and the nearest-neighbor pool is appropriately updated
    (as described above).</li>
    
    </ul>

    <p>This class includes functionality for both the Manager (MPI
    Rank == 0) and Worker (MPI Rank > 0) processes.  Necessary
    functionality to distinguish and operate either as Manager or
    Worker is already built into the class.  This class uses the
    MSTCache and MSTCluster classes to help in performing the various
    activities.  Refer to the documentation on the various methods for
    detailed description on their functionality and usage.</p>

*/
class NNMSTClusterMaker : public MSTClusterMaker {
    friend class ClusterMakerFactory;
public:
    /** The destructor.
        
        The destructor frees up any dynamic memory allocated by this
        object for its operations.
    */
    virtual ~NNMSTClusterMaker();

protected:
    /** A method to handle initialization tasks for the
		MSTClusterMaker.
		
        This method is called after the ESTs have been loaded and the
        clustering is just about to commence. This method calls the
        corresponding method in the base class to ensure base class
        initialization proceeds successfully (call to base class
        inititalizes the ESTAnalyzer and sets up estList pointer).  It
        then sets up the nearest-neighbor pool.

		\return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    virtual bool initialize();

    /** Helper method to compute \i local distance/similarity metrics
        and populate the information in the supplied smList.

        This method overrides the default implementation in the base
        class and is therefore called from the
        MSTClusterMaker::populateCache() method.  This method
        customizes the default operation by giving preference to the
        entries in the nearest-neighbor pool and searching for a match
        with those entries.  It performs an extensive comparison only
        if none of the entries in the nnPool yield a match.  The
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

    /**
       A simple method to print the current entries in the
       nearest-neighbor (NN) pool.
       
       This is utility method that can be used to dump the nnPool to a
       given output stream.  This is primarily used for
       troubleshooting this class.

       param[out] os The output stream to which the data is to be
       written.  Typically, this is just std::cout.
    */
    void printNNPool(std::ostream& os);
    
private:
    /** The constructor.
        
        The default constructor for this class.  The constructor is
        made private so that this class cannot be directly
        instantiated.  However, since the ClusterMakerFactory is a
        friend of this class, an object can be instantiated via the
        ClusterMakerFactory::create() method.

        \param[in,out] analyzer The EST analyzer to be used for
        obtaining similarity metrics between two ESTs.  This parameter
        is simply passed onto the base class.
    */
    NNMSTClusterMaker(ESTAnalyzer *analyzer);

	/** Helper method to perform a single analysis and appropriately
		add an entry to the smList.
	   
	*/
	bool computeSMEntry(const int estIdx,
                        const int otherEstIdx, SMList& smList,
						const float InvalidMetric,
						const bool addInvalidMetric,
                        const bool addToNNPool);

	/** Starting index of locally-owned ESTs that are managed by the
		nearest-neighbor pool (nnPool) in this class.

		This instance variable is a convenience instance variable that
		is used to track the starting index (zero-based) of the
		contiguous sub-set of ESTs that are operated on by this
		cluster maker.  Note that the set of ESTs to be processed are
		initially evenly divided between the parallel-processes used
		for clustering.  This value is set by the initialize method
		and is never changed during the life-time of this class. This
		value is used by the various methods in this class.
	*/
	int localESTStartIndex;

    /** A simple vector that indicates if a given EST index is present
        in the nearest-neighbor (nn) pool.

        This vector is setup in the initialize() method to contain a
        boolean flag for all locally-owned ESTs indicating if a given
        EST is in the pool. This pool is maintained locally on each
        parallel process (it does not contain information about all
        ESTs).
    */
	std::vector<bool> nnPool;
};

#endif
