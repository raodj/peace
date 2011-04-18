#ifndef FILTER_H
#define FILTER_H

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

#include "HashMap.h"
#include <vector>

// Forward declarations to keep compiler fast and happy
class RuntimeContext;
class ClusterMaker;
class ArgParser;
class EST;

/** \def FilteredESTList

    \brief Typedef for HashMap<int, std::vector<int> > FilteredESTList;

    This typedef is a convenience definition to define a hash map
    tracks the set of entries that were filetered out by this
    filter. The key to the hash map is the cluster ID to which
    filtered ESTs were added. The value in the hash map is a vector
    that contains the index of the ESTs that were added to the
    cluster.
*/
typedef HashMap<int, std::vector<int> > FilteredESTList;

/** The base class of all filters.
    
    <p>This class must be the base class of all filters in the
    system. This class provides some default functionality that can be
    readily used by the filters. This class enables the FilterChain to
    manage a list of filters and dispatch method calls to various
    filters.</p>

    <p>A set of filters (stored in the FilterChain class) are run on
    all entries in a FASTA file prior to commencement of the core
    clustering operation. The filters perform various validation
    operations to ensure ESTs are good prior to clustering. Filtering
    ensures that the overall quality of clustering provided by PEACE
    is good.</p>

    <p>Each filter object in the chain implements a specific type of
    filteration operation and ultimately returns an integer indicating
    the cluster to which an EST is to be assigned. If the cluster ID
    is -1, then that indicates that the EST must be subjected to
    regular clustering operations.</p>
*/
class Filter {
public:
    /** Add the set of command line parameters for this component.
        
        This method is invoked by the FilterChain component of PEACE
        (that contains this sub-component) when it is being requested
        for command line arguments.  This method is expected to add
        the various command line arguments that can be used to further
        customize the operations of just this filter.

        \note Derived classes must override this method.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser) = 0;
    
    /** Method to begin filter analysis (if any).
        
        This method is invoked (currently by the FilterChain) after
        command line arguments have been parsed but just before
        commencement of filtration.  First this method should validate
        and ensure it has the necessary command line parameters.
        Next, this method can load any additional information that may
        be necessary for a given filter.  In addition, it may perform
        any pre-processing as the case may be.

        \note Derived classes must override this method.
        
        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns with a
        non-zero error code.
    */
    virtual int initialize() = 0;

    /** Method to indicate completion of filter analysis.
        
        This method is invoked (currently by the FilterChain) after
        all the filteration operations have been successfully
        completed. This method typically performs any clean up
        operations that may be necessary.

        \note Derived classes must override this method.
    */
    virtual void finalize() = 0;

    /** Add cluster ID and indexes of ESTs filtered by this filter to
        a given list.

        This method is used to accumulate the set of ESTs that were
        filtered out by this filter into a single superList. The super
        list is then broadcasted to other processes for their
        reference.  The data from this filter is to be added to the
        superList in the following "flat" format (to ease broadcasting
        to other processes):

        <ol>

        <li>First the clusterID to which the ESTs were added must be
        appended to the superList.</li>

        <li>Next, the number of ESTs that were filtered out by this
        filter must be appended to the superList.</li>

        <li>Finally, the indexes of the ESTs that were filtered out by
        this filter must be appended to the superList.</li>
        
        </ol>

        \param[out] superList The vector to which the filter data is
        to be added. 
    */
    void addFilterData(std::vector<int>& superList) const;

    /** Helper method to process filter data from another process.

        This is a helper method to process a given list of entries
        (that was built by the addFilterData() method) obtained from
        another process. This list can contain entries from multiple
        filters (and not just one). This method method processes the
        list of entries and adds the indicated ESTs to corresponding
        clusters on local processes.  This ensures that filters that
        were independently applied on different parallel processes are
        consistently reflected on all processes participating in the
        clustering process.

        \note This method is present here so that the interface of
        this class is symmetric in the sense that the list generated
        by the addFilterData() method in this class is also processed
        by a complementary method in the same class.

        \note The method is static because it processes lists from
        multiple filters and not just one (so calling it on a single
        object does not make sense).  This is more of "purity" of API
        (if you are looking for some significance) rather than
        anything else.
        
        \param[in] superList The super list to be processed by this
        method.
        
		\param[in] clusterMaker The cluster maker object that is being
		used for analysis. This pointer can be NULL. If it is NULL,
		then fragments that were filtered out are not added to dummy
		clusters.
    */
    static void processFilterData(const std::vector<int>& superList,
                                  ClusterMaker *clusterMaker);
    
    /** Determine if the given EST passes this filter condition.
        
        This method can be used to determine if a given EST passes
        this filter and if it must be subjected to the core clustering
        operations.
        
        \param[in] est A pointer to the EST that must be subject
        to the filtering process.
        
        \return This method returns \c false if the EST must be
        subject to further filteration (or core EST analyis and
        clustering if this is the last filter in the chain).  If the
        specified EST is to be filtered out, then this method returns
        \c true.
    */
    bool applyFilter(const EST* est);

    /** Obtain human-readable name for this filter.

        This method must be used to obtain the human readable name set
        for this filter. This method essentially returns the value
        set when this filter class was instantiated.

        \return A human readable name associated with this filter.
    */
    const std::string& getName() const { return filterName; }
    
    /** The destructor.
        
        The destructor for the filter.  Currently it does not perform
        any special tasks.
    */
    virtual ~Filter();

    /** Method to display statistics regarding operation of this filter.

        This method can be used to obtain a dump of the statistics
        gathered regarding the operation of this filter.  The
        typical statistic generated by filters includes:

        <ul>

        <li>The number of times the filter was called.  More
        specifically this value indicates the number of times the \c
        applyFilter() method was invoked.</li>

        <li>The number of successful matches reported by this filter.
        This number indirectly indicates the number of times other
        filters were invoked.</li>
        
        </ul>

        \note Derived filter classes may override this method to
        display additional statistics. However, the additional
        information must be displayed after the base class method has
        completed its task.
        
        \param[out] os The output stream to which the statistics
        regarding the filter is to be dumped.
    */
    virtual void printStats(std::ostream& os) const;

    /** Method to obtain the count of times this filter was run.

        \return The number of times this filter was called.
    */
    inline int getRunCount() {
        return runCount;
    }

    /** Method to obtain the count of times this filter rejected (or
		filtered-out) an EST.

        \return The number of times calls to this filter rejected (or
        filtered-out) an EST.
    */
    inline int getFilterCount() {
        return filterCount;
    }
    
protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived Filter classes must be instantiated via the
        FilterFactory API methods.

        \param[in] filterName The human readable name for this filter.
        This name is used when generating errors, warnings, and other
        output messages for this filter.

        \param[in] runtimeContext The runtime context class from which
        the necessary context information (such as ESTAnalyzer to use
        or ClusterMaker to use) is to be extracted.  The pointer is
        saved so that the derived filter class can use it to perform
        specific operations.
    */
    Filter(const std::string& filterName, RuntimeContext *runtimeContext);

    /** Apply filter rules to determine if this EST should be filtered
        out.

        This method is invoked from the applyFilter() method to
        perform the actual filtering.  The filtering is performed
        on the given EST.

        \note Derived filter classes must override this method and
        provide a proper implementation.
        
        \param[in] est A pointer to an immutable EST object to which
        the filter rules are to be applied.

		\param[out] clusterID A non-negative integer value.  If this
        value is not -1, then it is used to place the EST into an
        artifically created cluster to help users identify such
        clusters.  The derived filter class is expected to
        appropriately fill in this parameter.

        \return This method returns \c false if the EST must be
        subject to further filteration (or core EST analyis and
        clustering if this is the last filter in the chain).  If the
        specified est is to be filtered out, then this method returns
        \c false.
    */
    virtual bool runFilter(const EST* est, int& clusterID) = 0;
    
    /** The name of this filter.
        
        This instance variable contains the human recognizable name
        for this filter.  This value is set when the filter is
        instantiated (in the constructor) and is never changed during
        the life time of this filter.  This information is used when
        generating errors, warnings, and other output messages.
    */
    const std::string filterName;

	/** Helper method to add a dummy cluster entry in the cluster maker.

		This method is a helper method that is used to add a dummy
		cluster entry to the ClusterMaker object that is being used
		for the current run of PEACE.  The ClusterMaker object is
		obtained from Filter::runtimeContext (that is set just before
		Filter::initialize() method is invoked by the FilterChain).
		This method uses the ClusterMaker::addDummyCluster() method to
		create the dummy cluster.
		
		\param[in] clusterInfo This the name or information about the
		cluster to be set when the cluster is created.

		\return This method returns the numerical ID of the dummy
		cluster.  This is a non-negative value (greater-than-or-equal
		to zero) if a dummy cluster was successfully created.
		However, if a valid ClusterMaker object was not present in the
		Filter::runtimeContext, then this method returns -1.
	*/
	int addDummyCluster(const std::string& clusterInfo);
	
    /** The runtime context set for this filter.

        This instance variable is initialized to refer to the
        RuntimeContext object passed to this filter for any operations
        that may be necessary. Note that this pointer cannot be change
        after it is set in the constructor.
    */
    RuntimeContext* const runtimeContext;
    
private:
    /** Variable to track the number of times this filter was run.

        This instance variable is used to track the number of times
        this filter was run. This variable is initialized to zero in
        the constructor.  It is incremented each time the
        applyFilter() method is invoked to run the the filter.
    */
    int runCount;

    /** Variable to track the number of times this filter filtered out
        an EST.

        This instance vairable tracks the number of times the filter
        filtered out an EST.  This value is incremented in the \c
        applyFilter() method each time the runFilter() method returns
        a \c non -1 value.
    */
    int filterCount;

    /** The list of ESTs filtered out by this filter.

        This instance variable is used to track the set of ESTs that
        were filtered out by this filter. A single filter can filter
        out ESTs and add them to different dummy
        clusters. Consequently a hash map is used to track the
        information.  The key to the hash map is the cluster ID to
        which filtered ESTs were added. The value in the hash map is a
        vector that contains the index of the ESTs that were added to
        the cluster. Entries to the hash map are added by the
        applyFilter method.
    */
    FilteredESTList filteredESTList;
    
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    Filter& operator=(const Filter& src);
};

#endif
