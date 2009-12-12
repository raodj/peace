#ifndef CLUSTER_MAKER_H
#define CLUSTER_MAKER_H

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

// Forward declaration to make compiler happy
class ESTAnalyzer;

/** The base class of all cluster makers.

    This class must be the base class of all cluster makers in the
    system. This class provides some default functionality that can be
    readily used by each cluster maker.
*/
class ClusterMaker {
public:
    /** Display valid command line arguments for this cluster maker.

        This method must be used to display all valid command line
        options that are supported by this cluster maker.  Note that
        derived classes may override this method to display additional
        command line options that are applicable to it.  This method
        is typically used in the main() method when displaying usage
        information.

        \note Derived cluster maker classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to display
        common options.
        
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

        \note Derived cluster maker classes <b>must</b> override this
        method to process any command line arguments that are custom
        to their operation.  When this method is overridden don't
        forget to call the corresponding base class implementation to
        display common options.
        
        \param[inout] argc The number of command line arguments to be
        processed.

        \param[inout] argc The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.  This method returns true if all arguments
        are consumed successfully.
    */
    virtual bool parseArguments(int& argc, char **argv);   
    
    /** Method to begin clustering.

        This method must be used to create clusters based on a given
        EST analysis method.  This method is a pure-virtual method.
        Therefore all cluster maker classes must override this method
        to perform all the necessary operations.
    */
    virtual int makeClusters() = 0;

    /** The destructor.

        The destructor frees memory allocated for holding any data in
        this base class.
    */
    virtual ~ClusterMaker();
    
protected:
    /** The default constructor.

        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived cluster maker classes must be instantiated via the
        ClusterMakerFactor API methods.

        \param[in] name The human readable name for this cluster
        maker.  This name is used when generating errors, warnings,
        and other output messages for this object.

        \param[inout] analyzer The EST analyzer to be used by this
        ClusterMaker for generating similarity metrics between two
        given ESTs.

        \param[in] refESTidx The reference EST's index in a given
        multi-FASTA file.  Index values start with 0 (zero).  The
        refESTidx is supplied as a global argument that is processed
        in the main() method.  This value is simply copied to the
        refESTidx member in this class.
        
        \param[in] outputFileName The file name to which output must
        be written.  If a valid output file is not specified, then
        results are written to standard output.  The outputFileName is
        simply copied to the outputFileName member object.
    */
    ClusterMaker(const std::string& name, ESTAnalyzer *analyzer,
                 const int refESTidx, const std::string& outputFileName);

    /** The name of this cluster maker.

        This instance variable contains the human recognizable name
        for this cluster maker.  This value is set when this object is
        instantiated (in the constructor) and is never changed during
        the life time of this object.  This information is used when
        generating errors, warnings, and other output messages.
    */
    const std::string name;

    /** The index of the reference EST in a given file. 

        This member object is used to hold the index of a reference
        EST in a given file.  The index values begin from 0 (zero).
        This member is initialized in the constructor and is never
        changed during the lifetime of this object.  This reference
        est index is typically the root node of any clustering
        operations that are performed.
    */
    const int refESTidx;

    /** The file to which results must be written.

        This member object is used to hold the file name to which all
        the clustering results are to be written.  This member is
        initialized to NULL.  However, the value is changed by the
        parseArguments method depending on the actual value specified
        by the user.
    */
    const std::string outputFileName;

    /** The analyzer to be used for generating EST similarity metrics.

        This pointer is used to hold a pointer to the EST analyzer
        that must be used for generating similarity metrics between
        two given pairs of ESTs.  This pointer is initialized when the
        object is instantiated and is never changed during the
        lifetime of this object.
    */
    ESTAnalyzer* const analyzer;
    
private:
    /** The set of common arguments for all cluster makers.

        This instance variable contains a static list of arguments
        that are common all the cluster makers.  The common argument
        list is statically defined and shared by all cluster maker
        instances.

        \note This makes cluster maker class hierarchy not MT-safe.
    */
    static arg_parser::arg_record commonArgsList[];

    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this
        object.
        
        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.
        
        \return Reference to this.
    */
    ClusterMaker& operator=(const ClusterMaker& src);
};

#endif
