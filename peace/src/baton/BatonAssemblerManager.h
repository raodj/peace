#ifndef BATON_ASSEMBLER_MANAGER_H
#define BATON_ASSEMBLER_MANAGER_H

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

#include "BatonAssembler.h"

// Some forward declarations to keep compiler fast and happy
class ContigMaker;

/** \file BatonAssemblerManager.h
    
    \brief The class definition for the manager (MPI rank 0) process
    of the parallel and distributed baton-based gene assembler.

    This file contains the class declaration for the manager process
    of the parallel baton assembler.  Additional documentation
    regarding the API provided by this class is available with the
    class and method documentation.
*/

/** The manager process in the parallel baton assembler.

    <p>This class serves as the manager process in the parallel and
    distributed architecture of the baton assembler.  There is a
    single, unique instance of this class (on the process with MPI
    rank zero) in the complete distributed system.  The manager is
    responsible for controlling and coordinating the processing
    performed by various distributed worker processes (if any).  If
    the baton assembler is being run as a single process, then the
    manager also performs the tasks that are typically designated to a
    worker process.</p>

    <p>This class is derived from the BatonAssembler class that serves
    as the common parent for both the worker and manager classes.  The
    shared base class implements the various common features of the
    worker and assembler.  This design streamlines the ability of the
    manager to also operate as a worker when the baton assembler is
    being run as a single process.</p>

    <p>This class overrides some of the core methods defined in the
    Assembler (grandparent) class, which is the primary API for all
    assemblers in PEACE, thereby completing the class hierarchy that
    implements all the API methods.  Consequently, this class is
    instantiable.  However, this class should not and cannot be
    directly instantiated.  Instead, the AssemblerFactory class must
    be used to create a suitable instance of this class.  The
    AssemblerFactory is responsible for ensuring that only a single
    instance of the BatonAssemblerManager class exists in all the
    processes associated with a given run of PEACE.</p>
*/
class BatonAssemblerManager : public BatonAssembler {
    friend class AssemblerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~BatonAssemblerManager();

    /** Primary method to perform cDNA assembly.
        
        This method is the primary interface method that is invoked
        (on each MPI process) to perform the complete assembly and
        output generation process.  This method is invoked once all
        the command line parameters are processed by the various
        subsystems constituting PEACE.  This method launches the
        various activities performed by the manager process.

		\note This method is invoked only after the initialize()
		method in the base class has successfully performed its
		operations.

        \return This method returns zero if the assembly process was
        successfully completed.  Otherwise this method returns a
        non-zero error code.
    */
    virtual int assemble();

	/** A method to handle initialization tasks.
		
		This method overrides the default implementation in the base
		class to initialize the progress file, if one has been
		specified.  This method is called after an assembler object
		has been instantiated but \i before the cDNAs (to be
		assembled) have been loaded.  This method performs the
		following tasks:

        <ol>

        <li> First it calls the base class implementation to loads the
		ESTs and other cDNA fragments to be processed from the
		specified input file.  The base class method initializes the
		baton analyzer.</li>

        <li>Next it attempts to open the progress file (if user has
        specified a file name) for writing progress information in the
        assemble() method.</li>
        
        </ol>
        
		\return This method returns zero on success. On errors, this
		method returns a non-zero value.
    */
    virtual int initialize();
	
protected:
    /** Helper method to initialize the set of ESTs to be assembled.

        This is a helper method that is invoked from the assemble
        method once to setup the list of ESTs to be assembled.  This
        method was introduced to streamline the code in the assemble
        method.  This list adds the index of all unprocessed (some
        ESTs are tagged as being processed by filters and other tools
        to indicate that they must be ignored) ESTs to the
        estsToProcess set.

        \note This method assumes that it will be called only once.
    */
    void setupESTsToBeProcessed();

	/** Helper method to add a list of alignment information to the
		contig maker and update pending set of fragments to be
		processed.

		This is a helper method that is called from the main
		assemble() method in this class after the manager has
		performed its part of checks.  This method is also called from
		the getWorkerAlignments() method when an alignment list is
		received from a worker.  For each alignment entry in the
		alignList parameter, this method adds the alignment to the
		contig maker. In addition, it removes the fragments added to
		the consensus from updates the estsToProcess set.

		\param[out] contig The contig maker object to which the
		alignment information in alignList must be added.

		\param[in] alignList The list containing alignment information
		to be added to the contig maker.  This parameter cannot be
		NULL.

		\param[in] listSize The number of elements in the alignList
		that must be added to the contig maker.
	*/
	void addToContig(ContigMaker& contig,
					 const AlignmentInfo* alingList,
					 const int listSize);

	/** Helper method to get alignment information from workers (if any).

		This method is a helper method that is called from the main
		align() method to obtain alignment information from the
		workers.  This methods performs various tasks only when:

		<ol>

		<li>MPI-based compilation has been enabled \b and </li>

		<li>There is at least one worker in the parallel run.</li>

		</ol>

		Each worker is expected to provide alignment information to
		the manager.  For each alignment list obtained from a worker,
		this method uses the addToContig() method to add the alignment
		list (obtained from a worker) to the overall contig maker.

		\note Every worker process always adds a dummy entry to the
		alignment list at the end.  This last entry is ignored.
		However, having this dummy entry guarantees that every worker
		will always deliver an alignment list to the manager,
		streamlining the logic of handling distributed workers.
		
		\param[out] contig The contig maker to which the alignment
		information obtained from a worker is to be added.
	*/
    void getWorkerAlignments(ContigMaker& contig);

    /** Helper method to write progress information to a given data file.

        This method is a helper method that is invoked from the main
        assemble method to report progress logs as ESTs are analyzed
        and assembled.  This method cuts logs only if the progressFile
        is valid.  It logs the number of ESTs processed and the total
        number of fragments using the values returned by
        EST::getESTCount() and the number of entries left in
        estsToProcess.
    */
    void reportProgress();
    
private:
    /* The constructor for this class.
       
       The constructor has been made private so that this class cannot
       be directly instantiated.  However, notice that the
       AssemblerFactory is a friend of this class; therefore it can
       instantiate a BatonAssemblerManager object.  Accordingly, the
       AssemblerFactory::create() method must be used to instantiate
       an object of this class.

       \param[in] outputFile The name of the output file to which the
       assembled information is to be written.  This value is simply
       passed onto the base class as the output file name.  Refer to
       base class documentation for details on usage of this
       parameter.
    */
    BatonAssemblerManager(const std::string& outputFileName);

    /** The set of cDNA fragments to be processed by any process.

        This vector is used to maintain the index of the cDNA
        fragments to be assembled across all MPI processes.  This list
        is initialized first-thing in the assemble method.  It is
        filled-in with the index of the complete set of cDNA fragments
        to be assembled.  As cDNA fragments are assembled the
        corresponding index values removed from this list by the
        removeProcessedEntries() method.  When this list becomes
        empty, that means all the fragments have been assembled and it
        is time to windup assembly.
    */
    std::set<int> estsToProcess;

    /** File stream to log progress information.
        
        This output stream is created in the initialize() method if
        the progFileName is specified.  The progress information is
        generated by the updateProgress method.
    */
    std::ofstream progressFile;
};

#endif
