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
#include "BatonListCache.h"
#include "DefaultSequenceAligner.h"
#include "MPIHelper.h"

#include <set>
#include <fstream>

// Forward declarations
class BatonAnalyzer;
class ContigMaker;

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

    /** Wind-up the operations of this component.

        This method is invoked when the SubSystem (that logically owns
        this component) is being finalized.  This method currently
        prints MPI statistics.
    */
    virtual void finalize();

    
    /** Common method to determine the cDNA fragment to be used as the
        initial reference for contig generation.

        This is a common method that is invoked by the manager and
        worker(s) to determine the reference fragment to be used as
        the initial seed for contig creation.  All the parallel
        processes collaboratively compute the next available cDNA
        fragment to be used as the initial seed.

        \note Currently, this method selects the longest, un-processed
        cDNA fragment as the candidate.
    */
    int getReferenceEST();
    
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

    /** Compute alignment between a reference EST and the subset of
        ESTs owned by this assembler process.

        This method is a common helper method that is shared by the
        manager and worker processes.  In this context the term
        'local' is used to refer to the local process on which this
        method is invoked.  This method is repeatedly invoked (with
        different reference fragments) to compute alignment between a
        given reference sequence (that has already been set/computed
        by the ContigMaker) and the subset of entries owned by this
        assembler process.
        
        \note The subset of fragments owned by this process is
        computed via a call to the getOwnedESTidx() method.

        \note Fragments that were aligned are flagged as having been
        processed by this method.
        
        \param[out] contigMaker The contig maker object to be used by
        this method to form the contig.

		\param[in] isManager This parameter must be set to true if
		this method is being called by the manager process. The
		workers must set this flag to false.

		\return This method returns the number of ESTs that were
		processed in this call across all parallel processes. This
		information is typically used by the manager to report
		progress information.
    */
    int localAssembly(ContigMaker& contigMaker, const bool isManager);

    /** Helper method to create a new Contig class and add it to the
        list of contigs created by this assembler.

        This is a common (used by both manager and worker) to create
        and add a new contig. The contig object is created using the
        information from supplied ContigMaker (parameter). The created
        contig is added to the contigList member (defined in the base
        class).  This method is created from the localAssembly()
        method once the contig is ready.

        \note All parallel processes must call this method to ensure
        that contigs are consistently added all processes to preserve
        all the necessary information about the contig.

        \param[in] contigMaker The contig maker object (that has
        finished building a contig) from which contig information is
        to be obtained to create a new contig.

        \param[in] isManager If this flag is \c true, then this method
        also populates the consensus sequence and nucleotide
        distribution data into the new contig created by this method.
    */
    void createAddContig(const ContigMaker& contigMaker, const bool isManager);

    /** Helper method to obtain total number of ESTs processed <i>only
        at the manager process</li> by all parallel-processes.
        
        This method is a helper method that can be used to obtain the
        total number of ESTs have been marked has having been
        processed.  This method is primarily used to collate data from
        all the parallel processes involved in assembly.
        Consequently, this method must be invoked in a
        coordinated-manner on all the parallel-processes. This method
        adds up the number of processed ESTs at each process and
        returns the total at the given process ID. At all other
        processes the local count is returned.

        \return The sum of number of ESTs that have already been
        processed at the process with MPI-rank == procID.  At all
        other processes this method returns the local count of the
        number of ESTs processed.
    */
    int getESTsProcessed(const int procID = MANAGER_RANK) const;
    
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
    /** Helper method to grow a given contig by matching unprocessed
        ESTs with the contig's root cDNA fragment.

        This method is a helper method that is invoked (one or more
        times) from the localAssembly() method in this class. This
        method uses batons to identify batons with sufficient
        similarity and then constructs segments representing regions
        of overlap. The segments are added to the supplied contig
        maker.
        
        \param[out] contig The contig maker object to be used to
        explore and build the contig. This contig can either be the
        primary starting contig or a contig to its left or right.

        \return This method returns the number of cDNA fragments that
        were added to the contig maker.  If no ESTs were added to the
        contig maker, then this method returns zero.        
    */
    int growContig(ContigMaker& contig);
};

#endif
