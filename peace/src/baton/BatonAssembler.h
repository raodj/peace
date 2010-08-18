#ifndef BATON_ASSEMBLER_H
#define BATON_ASSEMBLER_H

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

#include "Assembler.h"
#include "AlignmentInfo.h"

#include <set>
#include <fstream>

/** \file BatonAssembler.h
    
    \brief The common base class for the parallel and distributed
    baton-based gene assembler.

    This file contains the common base class declaration for the baton
    assembler.  Additional documentation regarding the API provided by
    this class is available with the class and method documentation.
*/

/** A common base class for the baton assembler.

	<p>The (at least conceptual) baton assembler is a parallel and
    distributed algorithm whose functionality can be broadly
    classified into two main categories, namely: a manager process and
    a set of worker processes.  Accordingly, the baton assembler has
    been implemented using two different classes, namely: BatonManager
    and BatonWorker.  This class provides a common base class for
    housing functionality that is common to the aforementioned two
    child classes.</p>
	
    <p>This assembler uses the concept of batons to perform gene
    assembly.  The baton approach essentially uses the maximum number
    of identical batons between two given cDNA fragments to determine
    if the pair are related and can be assembled together.</p>

    <p>In addition to housing common functionality, this class also
    implements the interfaces required by the more general purpose
    Assembler base class.  However, this class cannot be directly
    instantiated.  Instead one of the derived child classes must be
    instantiated and used.</p>
*/
class BatonAssembler : public Assembler {
public:
    /** The set of tags exchanged between manager and worker processes
        
        This enum provides meaningful names to the various tags
        (integers) exchanged between the manager and worker processes
        participating in baton-based assembly in a
        parallel/distributed manner.
    */
    enum MessageTags{UNKNOWN_REQUEST, ALIGN_INFO_LIST};
  
    /** The destructor.
        
        The destructor frees up any dynamic memory allocated by this
        class.  The destructor also deletes the instance of the
        BatonAnalyzer class created in the constructor.
    */
    virtual ~BatonAssembler();
    
    /** Display valid command line arguments for this assembler.
	
        This is used to display all valid command line options that
        are supported by this assembler.  Many of the options for
        specifying window size, <i>n</i>-mer size, threshold,
        permitted errors, and sufficiently good alignment score is
        delegated to the BatonAnalyzer member in this class.

        \note The Assembler base class requires that derived EST
        assembler classes <b>must</b> override this method to display
        help for their custom command line arguments.  When this
        method is overridden don't forget to call the corresponding
        base class implementation to display common options.
        
        \param[out] os The output stream to which the valid command
        line arguments must be written.
    */
    virtual void showArguments(std::ostream& os);

    /** Process command line arguments.

        This method is used to process command line arguments specific
        to this assembler.  This method is typically used from the
        core system just after an derived assembler has been
        instantiated (MPI (if used) has already been initialized).
        This method must consumes all valid command line arguments
        applicable to the implementation.  If the command line
        arguments were valid and successfully processed, then this
        method returns \c true.

        \note Derived assembler classes are expected to override this
        method to process any command line arguments that are custom
        to their operation.  When this method is overridden don't
        forget to call the corresponding base class implementation to
        display common options.
        
        \param[in,out] argc The number of command line arguments to be
        processed.

        \param[in,out] argv The array of command line arguments.

        \return This method returns \c true if the command line
        arguments were successfully processed.  Otherwise this method
        returns \c false.  This method returns true if all arguments
        are consumed successfully.
    */
    virtual bool parseArguments(int& argc, char **argv);

	/** A method to handle initialization tasks.
		
		This method overrides the default implementation in the base
		class to initialize the baton analyzer.  This method is called
		after an assembler object has been instantiated but \e before
		the cDNAs (to be assembled) have been loaded.  This method
		performs the following tasks:

        <ol>

        <li> First it calls the base class implementation to loads the
		ESTs and other cDNA fragments to be processed from the
		specified input file.</li>

        <li>Next it initializes the BatonAnalyzer (via method that
        does not load data).</li>
        
        </ol>
        
		\return This method returns zero on success. On errors, this
		method returns a non-zero value.
    */
    virtual int initialize();
		
protected: 
    /* The constructor for this class.
       
       The constructor is made protected so that this class cannot be
       directly instantiated.  Instead one of the derived classes must
       be suitably instantiated.

       \param[in] outputFile The name of the output file to which the
       assembly information is to be written.  If this parameter is
       the empty string then output is written to standard output.
       This value is simply passed onto the base class.
    */
    BatonAssembler(const std::string& outputFileName);

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

		\note Currently, this method has exact implementation as in
		MSTClusterMaker::getOwnedESTidx. This method has been
		copy-pasted so that filters can operate on their own different
		sub-set of ESTs if they choose. Maybe the method can be
		combined together.
    */
    static void getOwnedESTidx(int& startIndex, int& endIndex);

    /** The baton analyzer that does pair-wise alignment.

        This instance variable provides a reference to the baton
        analyzer that performs the core task of analyzing a pair of
        cDNA sequences and providing a suitable alignment for them,
        assuming that the two sequences are sufficiently similar to
        warrant alignment.  This pointer is initialized in the
        constructor (via the ESTAnalyzerFactory) and deleted in the
        destructor.
    */
    BatonAnalyzer* const analyzer;

    /** Compute alignment between refESTidx and the subset of ESTs
        owned by this assembler process.

        This method is a common helper method that is shared by the
        manager and worker processes.  This method is repeatedly
        invoked (with different reference fragments) to compute
        alignment between a given reference sequence and the subset of
        entries owned by this assembler process.
        
        \note The subset of fragments owned by this process is
        computed via a call to the getOwnedESTidx method in this
        class.

        \note Fragments that were aligned are flagged as having been
        processed by this method.
        
        \param[in] refESTidx The index of the reference sequence to be
        used for exploring and discovering similar fragments that can
        be aligned together.

        \param[out] alignList This vector is populated with the list
        of alignments that were discovered in the subset of fragments
        owned by this process.  Note that this vector is not cleared
        by this method.  New alignment information is simply added to
        this list.
    */
    void localAssembly(const int refESTidx, AlignmentInfoList& alignList);

    /** Name of file to report progress in during assembly.

        This command line argument provides the name of the log file
		where progress information is to be written. The progress
		information is in the form: \#estsProcessed, \#ests. The file
		name is specified via a command line argument \c --progress.
		This value is used only by the BatonAssemblerManager class.
		Possibly this parameter must be moved down to the child class.
    */
    static char *progFileName;
	
private:
    /** The set of arguments specific to the baton assembler.

        This instance variable contains a static list of command line
        arguments that are specific only to the baton assembler class.
        This argument list is statically defined and shared by all
        instances of this class.

        \note Use of static arguments and parameters makes this class
        hierarchy not MT-safe.
    */
    static arg_parser::arg_record argsList[];
};

#endif
