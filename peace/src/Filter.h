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

#include "arg_parser.h"
#include "HashMap.h"

#include <vector>

// Forward declarations to keep compiler fast and happy
class ClusterMaker;

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
    /** Display valid command line arguments for this filter.
        
        This method must be used to display all valid command line
        options that are supported by this filter.  Note that derived
        classes may override this method to display additional command
        line options that are applicable to it.  This method is
        typically used in the main() method when displaying usage
        information.

        \note Derived filter classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to display
        common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os) = 0;

    /** Process command line arguments.
        
        This method is used to process command line arguments specific
        to this filter.  This method is typically used from the main
        method just after the filter has been instantiated.  This
        method consumes all valid command line arguments.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.

        \note Derived filter classes <b>must</b> override this method
        to process any command line arguments that are custom to their
        operation.  When this method is overridden don't forget to
        call the corresponding base class implementation to display
        common options.
        
        \param[in,out] argc The number of command line arguments to be
        processed. This value is updated when valid command line
        arguments are consumed by the filter. 
        
        \param[in,out] argv The array of command line arguments. The
        number of entries in this array are modified and updated when
        valid arguments are consumed by the filter.
        
        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.
    */
    virtual bool parseArguments(int& argc, char **argv) = 0;
    
    /** Method to begin filter analysis (if any).
        
        This method is invoked just before commencement of filtration.
        This method typically loads additional information that may be
        necessary for a given filter.  In addition, it may perform any
        pre-processing as the case may be.

        \note Derived classes must override this method.
        
        \return If the initialization process was sucessful, then this
        method returns 0.  Otherwise this method returns with a
        non-zero error code.
    */
    virtual int initialize() = 0;

    /** Method to indicate completion of filter analysis.
        
        This method is invoked after all the filteration operations
        have been successfully completed. This method typically
        performs any clean up operations that may be necessary.

        \note Derived classes must override this method.
    */
    virtual void finalize() = 0;

    /** Add cluster ID and indexes of ESTs filtered by this filter to
        a list.

        This method is used to accumulate the set of ESTs that were
        filtered out by this filter into a single superList. The super
        list is then broadcasted to other processes for their
        reference.  The data from this filter is to be added to the
        superList in the following format (to ease broadcasting to
        other processes):

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
        filters. This method method processes the list of entries and
        adds the indicated ESTs to corresponding clusters on local
        processes.  This ensures that filters that were independently
        applied on different parallel processes are consistently
        reflected on all processes participating in the clustering
        process.

        \note This method is present here so that the interface of
        this class is symmeetric in the sense that the list generated
        by the addFilterData method is also processed by this method.

        \note The method is also static because it processes lists
        from multiple filters and not just one.  This is more of
        "purity" of API and its significance rather than anything
        else.

        \param[in] superList The super list to be processed by this
        method.

		\param[in] clusterMaker The cluster maker object that is being
		used for analysis.
    */
    static void processFilterData(const std::vector<int>& superList,
                                  ClusterMaker *clusterMaker);
    
    /** Determine if the given EST passes this filter condition.
        
        This method can be used to determine if a given EST passes
        this filter and if it must be subjected to the core clustering
        operations.
        
        \param[in] otherEST The index (zero based) of the EST that
        must be subject to which the reference EST is to be compared.
        
        \return This method returns -1 if the est must be subject to
        further filteration or core EST analyis and clustering.  If
        the specified est is to be filtered out, then this method
        returns a non-zero integer value. This value is used to place
        the EST into an artifically created cluster to help users
        identify such clusters.
    */
    int applyFilter(const int otherEST);

    /** Obtain human-readable name for this filter.

        This method must be used to obtain the human readable name set
        for this filter. This method essentially returns the value
        set when this filter class was instantiated.

        \return A human readable name associated with this filter.
    */
    const std::string& getName() const { return filterName; }
    
    /** The destructor.
        
        The destructor for the filter.
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

        \param[in] clusterMaker The cluster maker class that is being
        used for analysis.  This class can be used by this filter for
        performing specific operations.
    */
    Filter(const std::string& filterName, ClusterMaker *clusterMaker);

    /** Apply filter rules to determine if this EST should be filtered
        out.

        This method is invoked from the applyFilter() method to
        perform the actual filtering.  The filtering is performed
        on the given EST.

        \note Derived filter classes must override this method and
        provide a proper implementation.
        
        \param[in] estIndex The index (zero based) of the EST to which
        the filter rules are to be applied.

        \return This method returns -1 if the est must be subject to
        further filteration or core EST analyis and clustering.  If
        the specified est is to be filtered out, then this method
        returns a non-zero integer value. This value is used to place
        the EST into an artifically created cluster to help users
        identify such clusters.
    */
    virtual int runFilter(const int estIndex) = 0;
    
    /** The name of this filter.
        
        This instance variable contains the human recognizable name
        for this filter.  This value is set when the filter is
        instantiated (in the constructor) and is never changed during
        the life time of this filter.  This information is used when
        generating errors, warnings, and other output messages.
    */
    const std::string filterName;

    /** The cluster maker set for this filter.

        This instance variable is initialized to refer to the
        top-level cluster maker class to be used by this filter for
        any operations that may be necessary. Note that this pointer
        cannot be change after it is set in the constructor.
    */
    ClusterMaker* const clusterMaker;
    
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
