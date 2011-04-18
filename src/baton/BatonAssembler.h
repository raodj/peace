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
#include "BatonListCache.h"
#include "DefaultSequenceAligner.h"

#include <set>
#include <fstream>

// Forward declarations
class BatonAnalyzer;

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
        SequenceAnalyzer object created in the initialize() method.
    */
    virtual ~BatonAssembler();

	/** Add valid command line arguments for this assembler.
		
        This method must be used to add all valid command line options
        that are supported by this assembler.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.
		
        \note Derived assembler classes <b>must</b> override this
        method to display help for their custom command line
        arguments.  When this method is overridden don't forget to
        call the corresponding base class implementation to add common
        options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

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
        
		\return This method returns \true on success. On errors, this
		method returns \c false.
    */
    virtual bool initialize();
		
protected: 
    /* The constructor for this class.
       
       The constructor is made protected so that this class cannot be
       directly instantiated.  Instead one of the derived classes must
       be suitably instantiated.
    */
    BatonAssembler();

    /** The sequence aligner that does pair-wise alignment.

        This instance variable provides a reference to the sequence
        analyzer that performs the core task of analyzing a pair of
        cDNA sequences and providing a suitable alignment for them,
        assuming that the two sequences are sufficiently similar to
        warrant alignment.  
    */
    DefaultSequenceAligner aligner;

    /** Compute alignment between a given EST and the subset of ESTs
        owned by this assembler process.

        This method is a common helper method that is shared by the
        manager and worker processes.  In this context the term
        'local' is used to refer to the local process on which this
        method is invoked.  This method is repeatedly invoked (with
        different reference fragments) to compute alignment between a
        given reference sequence and the subset of entries owned by
        this assembler process.
        
        \note The subset of fragments owned by this process is
        computed via a call to the getOwnedESTidx() method.

        \note Fragments that were aligned are flagged as having been
        processed by this method.
        
        \param[in] refEST A pointer to the reference EST to be used
        for exploring and discovering similar fragments that can be
        aligned together.

        \param[out] alignList This vector is populated with the list
        of alignments that were discovered in the subset of fragments
        owned by this process.  Note that this vector is not cleared
        by this method.  New alignment information is simply added to
        this list.
    */
    void localAssembly(const EST* refEST, AlignmentInfoList& alignList);

	/** A shared BatonList cache to cache baton lists.

		This object is used to contain a shared list of BatonList
		objects that are cached and reused during various analysis.
		Entries in the cache are created on demand.  In addition to
		the BatonAssembler hierarchy, the BatonAnalyzer also uses this
		cache for analysis.  The BatonListCache adds a couple of
		command line arguments to customize its operations.
	*/
	BatonListCache blCache;
	
private:
    // Currently this class does not have any private members
};

#endif
