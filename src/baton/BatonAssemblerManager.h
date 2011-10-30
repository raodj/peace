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
		has been instantiated but \e before the cDNAs (to be
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
        
		\return This method returns \c true on success. On errors, this
		method returns \c false.
    */
    virtual bool initialize();
	
protected:
	// Currently this class does not have any protected members
	
private:
    /* The constructor for this class.
       
       The constructor has been made private so that this class cannot
       be directly instantiated.  However, notice that the
       AssemblerFactory is a friend of this class; therefore it can
       instantiate a BatonAssemblerManager object.  Accordingly, the
       AssemblerFactory::create() method must be used to instantiate
       an object of this class.
    */
    BatonAssemblerManager();

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
};

#endif
