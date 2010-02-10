#ifndef LENGTH_FILTER_H
#define LENGTH_FILTER_H

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

/** A simple filter to filter out reads shorter than a given length.
    
	<p>This class provides a simple filter that can be used to filter
	out ESTs that are shorter than a given number of nucleotides.
	This filter can be enabled by specifying its name, namely \c
	lengthFilter, in the filter chain (that is specified as the
	command line argument to PEACE). If this filter is specified then
	the lengthFilterSize command line argument is used to determine
	the threshold length value (ESTs shorter than this value will be
	filtered out).</p>

    <p>This filter creates an dummy cluster with the meta name "Short
    ESTs (filtered by LengthFilter)" and adds any ESTs filtered out by
    this EST to that cluster.</p>

    <p>This class has been developed by extending the Filter base
    class and implementing the necessary API methods specified by the
    base class.  This enables the LengthFilter to be used in the
    FilterChain along with other filters to filter out ESTs. In
    addition, note that this class cannot be directly
    instantiated. Instead, the FilterFactory::create() method must be
    used to obtain an instance of this class.</p>
*/
class LengthFilter : public Filter {
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
        This method adds a new dummy cluster using the ClusterMaker
        API methods.  The dummy cluster is used to classify ESTs that
        were filtered out by this filter.

        \return This method returns zero to indicate that
        initialization was completed successfully. On errors (that is
        a cluster could not be created) then this method returns a
        non-zero error code.
    */
    virtual int initialize();

    /** Method to indicate completion of filter analysis.
        
        This method is invoked after all the filteration operations
        have been successfully completed.  Currently method this
        method does not have any specific tasks to perform.
    */
    virtual void finalize() {}
    
    /** The destructor.
        
        The destructor for the filter. The destructor currently has no
        specific tasks to perform as this filter does not use any
        dynamic memory.
    */
    virtual ~LengthFilter() {}
    
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
    LengthFilter(ClusterMaker *clusterMaker);
    
    /** Apply filter rules to determine if this EST should be filtered
        out.

        This method is invoked from the applyFilter() method to
        perform the actual filtering.  The filtering is performed on
        the given EST in the following manner:

        <ol>

        <li>This method first obtains the specified EST from the list
        of ESTs.</li>

        <li>If the nucleotide sequence of this EST is shorter than \c
        minESTLength then this method filters out this EST and returns
        the ID of the cluster to which this EST must be added. The ID
        of the cluster is determined in the LengthFilter::initialize()
        method.</li>

        \return This method returns -1 if the est must be subject to
        further filteration or core EST analyis and clustering.  If
        the specified est is to be filtered out, then this method
        returns a non-zero integer value. This value is used to place
        the EST into an artifically created cluster to help users
        identify such clusters.
    */
    virtual int runFilter(const int otherEST);
    
private:
    /** Variable to track the ID of the cluster to which short ESTs
        must be added.

        This instance variable is used to track the ID of the cluster
        to which short ESTs must be added by this filter. This value
        is set in the LengthFilter::initialize() method.
    */
    int clusterID;

    /** The set of arguments specific to this filter.

        This instance variable contains a static list of arguments
        that are specific only to this filter class.  This argument
        list is statically defined and shared by all instances of this
        class.

        \note Use of static arguments and parameters renders this
        filter class not to be MT-safe.
    */
    static arg_parser::arg_record argsList[];

    /** The minmum length threshold for ESTs to be used by this
        filter.

        This instance variable contains the length of the shortest EST
        that must be permitted through by this filter. All ESTs
        shorter than this length will be filtered out. The default
        length is 50. However, this value can be overridden by using
        the \c --minESTLen command line argument.
    */
    static int minESTLen;
    
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    LengthFilter& operator=(const LengthFilter& src);
};

#endif
