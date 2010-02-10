#ifndef LC_FILTER_H
#define LC_FILTER_H

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
#include "Filter.h"

// Forward declarations to keep compiler fast and happy
class ClusterMaker;

/** \typedef std::pair<int, SMList> DummyESTInfo

    \brief Shortcut typedef for std::pair<int, int>

    This typedef is a shortcut for referring to a pair of integers
    that hold information about dummy ESTs created, used, and finally
    removed by this filter. The pair of information is used as
    follows:

    <ol>

    <li>The first entry is an \c int that indicates the index of the
    dummy EST in the global EST list. </li>

    <li>The second entry is another integer hat contains the ID of the
    cluster to which other ESTs that are sufficiently close to the
    dummy EST must be added.</li>

    </ol>
*/
typedef std::pair<int, int> DummyESTInfo;

/** A filter to weed out reads with Low Complexity (LC) sections.
    
	<p>This class provides a filter that can be used to filter out
	ESTs that contain regions of Low Complexity (LC) reads in them.
	This filter is needed When clustering FASTA data (that contains
	low complexity reads) because the LC sections provide a "false"
	relationships between ESTs giving raise to very large clusters.
	These large clusters are created because transitive relationships
	are established between ESTs due to low complexity reads.</p>

	<p>In order to avoid super-clusters that get formed due to low
	complexity reads, PEACE adds dummy ESTs (such as: one with all \c
	"AAAAA...." and another with all \c "CCCCCC...") based on the
	pattern ESTs specified by the user. The length of the dummy ESTs
	are twice the length of the largest window used for analysis. The
	ESTs are subjected to the same analysis and ESTs that are
	sufficiently similar to the dummy entries are filtered out.</p>
	
    <p>This filter creates several dummy clusters (one per pattern
    sepcified) with the meta name "Low Complexity ESTs (filtered by
    LCFilter Pattern AA)" and adds any ESTs filtered out by this EST
    to the appropriate cluster.</p>

    <p>This class has been developed by extending the Filter base
    class and implementing the necessary API methods specified by the
    base class.  This enables the LCFilter to be used in the
    FilterChain along with other filters to filter out ESTs. In
    addition, note that this class cannot be directly
    instantiated. Instead, the FilterFactory::create() method must be
    used to obtain an instance of this class.</p>
*/
class LCFilter : public Filter {
    friend class FilterFactory;
public:
    /** Display valid command line arguments for this filter.
        
        This method must be used to display all valid command line
        options that are supported by this filter.  This method
        overrides the corresponding method in the base class API. This
        method is typically used in the main() method when displaying
        usage information.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);
    
    /** Process command line arguments.
        
        This method is used to process command line arguments specific
        to this filter.  This method is typically used from the main
        method just after the filter has been instantiated.  This
        method consumes all valid command line arguments.  If the
        command line arguments were valid and successfully processed,
        then this method returns \c true.
        
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
    virtual bool parseArguments(int& argc, char **argv);

    /** Method to begin filter analysis (if any).
        
        This method is invoked just before commencement of filtration.
        This method creates a set of dummy ESTs and corresponding
        dummy clusters for filtering out entries with low complexity
        regions. It uses the patternList to add dummy ESTs.

        \return This method returns zero to indicate that
        initialization was completed successfully. On errors (that is
        a dummy EST could not be created) then this method returns a
        non-zero error code.
    */
    virtual int initialize();

    /** Method to indicate completion of filter analysis.
        
        This method is invoked after all the filteration operations
        have been successfully completed.  This method removes all the
        dummy ESTs that were created and added by this
        filter.
    */
    virtual void finalize();
    
    /** The destructor.
        
        The destructor for the filter. The destructor currently has no
        specific tasks to perform as this filter does not use any
        dynamic memory.
    */
    virtual ~LCFilter() {}
    
protected:
    /** The default constructor.
        
        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead an instance
        should be created via a suitable call to the FilterFactory API
        method(s).
        
        \param[in] clusterMaker The cluster maker class that is being
        used for analysis.  This parameter is simply passed onto the
        base class for its use. It is used the initialize method to
        create a dummy cluster for use by this filter.
    */
    LCFilter(ClusterMaker *clusterMaker);
    
    /** Apply filter rules to determine if a given EST should be
        filtered out.

        This method is invoked from the applyFilter() method to
        perform the actual filtering.  The filtering is performed on
        the given EST in the following manner:

        <ol>

        <li>This method obtains the analyzer from the cluster maker
        and sets the given estIdx as the reference EST for
        analysis.</li>

        <li>For each dummy EST in \c dummyESTList this method performs
        the following tasks:

        <ol>

        <li>It uses the analyzer (in the \c clusterMaker) to compute a
        metric between the dummy EST and the given \c estIdx.</li>

        <li>If the metric indicates that the two ESTs are closely
        related (based on the threshold specified by the user) then
        the EST is filtered out and added to the corresponding dummy
        cluster.</li>

        </ol>

        </li>

        \param[in] estIdx The index of the EST to be tested and
        filtered by this method.
        
        \return This method returns -1 if the est must be subject to
        further filteration or core EST analyis and clustering.  If
        the specified est is to be filtered out, then this method
        returns a non-zero integer value. This value is used to place
        the EST into an artifically created cluster to help users
        identify such clusters.
    */
    virtual int runFilter(const int estIdx);
    
    /** Add a dummy entry with given sequence and length to the list
        of ESTs.
        
        This method is invoked from the initialize method to add a
        dummy EST. Typically two dummy ESTs (one all \c "AAAAA..."
        and another one with all \c "CCCC....") are added. This is 

        \note This method approriately adds an entry to the \c
        dummyESTList.
        
        \param[in] fastaID The fasta ID to be assigned to the dummy
        EST. This FASTA ID is not really useful is included for
        completeness.

        \param[in] seq The nucleotide sequences to be repeated in
        order to generate the complete FASTA sequence for the dummy
        EST. This sequence must be at least one character in
        length. Only valid nucleotide base pairs must be 

        \param[in] length The minimum length of the generated
        sequence. Note that if the length and the pattern in the
        squence are not integral multiples of each other, then this
        method may generate slightly longer ESTs.
    */
    virtual void addDummyEntry(const std::string& fastaID,
                               const std::string& seq,
                               const int length);

private:
    /** A list containing information about dummy ESTs.

        This list is used to hold the core information about dummy
        ESTs created, used, and (finally) removed by this filter.  For
        each pattern specified by the user (as a command line
        argument), a dummy EST is added (to the global list of ESTs)
        and this list maintains information about them. Entries are
        added to this list in the initialize() method. The entries are
        used on the runFilter() method. The finalize() method clears
        out the dummy ESTs and the entries in this list.
    */
    std::vector<DummyESTInfo> dummyESTList;
    
    /** The set of arguments specific to this filter.

        This instance variable contains a static list of arguments
        that are specific only to this filter class.  This argument
        list is statically defined and shared by all instances of this
        class.

        \note Use of static arguments and parameters renders this
        filter class not to be MT-safe.
    */
    static arg_parser::arg_record argsList[];
    
    /** The list of patterns to be used for generating dummy ESTs.
        
        This variable is used to refer to the set of patterns that
        must be used to generate the dummy ESTs for identifying low
        complexity regions. The list is in the form "A,C,AG,TC" --
        that is it contains a comma separated list of patterns. The
        patterns are repeated to generate dummy entries.  The pattern
        can eb set via the \c --lcPatterns command line argument.
    */
    static char* patternList;

    /** The default pattern used by this filter.

        This constant defines the default pattern that is used by this
        filter to generate dummy ESTs for identifying low complexity
        regions.  The current default pattern is "A,C". This causes
        two dummy ESTs (one with all \c AAAAA... and another EST with
        all \c CCCC...) to be created and used for filtering.
    */
    static char DefaultPatternList[];
    
    /** The filter's similarity/distance threshold metric.

        This variable is used to contain the similarity/distance
        metric to be used as the threshold value. This value is
        compared with the metric provided by the ESTAnalyzer to
        determine if a given EST contains a low complexity
        section. This value is an important metric. This value is set
        via the \c --lcThreshold command line argument.
    */
    static int threshold;
    
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    LCFilter& operator=(const LCFilter& src);
};

#endif
