#ifndef BATON_ASSEMBLER_WORKER_H
#define BATON_ASSEMBLER_WORKER_H

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

/** \file BatonAssemblerWorker.h
    
    \brief The class definition for the worker (MPI rank &gt 0)
    process of the parallel and distributed baton-based gene
    assembler.

    This file contains the class declaration for the worker process
    of the parallel baton assembler.  Additional documentation
    regarding the API provided by this class is available with the
    class and method documentation.
*/

/** The worker process in the parallel baton assembler.

    <p>This class serves as the worker process in the parallel and
    distributed architecture of the baton assembler.  There are
    multiple worker processes (however, each process has only a single
    instance of this class) in a parallel run of the baton assembler.
    The workers are run on MPI processes that have a non-zero MPI
    rank.  The worker is responsible for coordinating and assisting
    the manager in baton-based analysis and assembly.</p>

	\note If the baton assembler is being run as a single process,
    then the manager performs the tasks that are typically designated
    to a worker process.  Consequently, in a single process scenario
    this worker class is not instantiated.

    <p>This class is derived from the BatonAssembler class that serves
    as the common parent for both the manager and worker classes.  The
    shared base class implements the various common features of the
    worker and manager.  Consequently, many of the methods in the
    worker simply delegate tasks to various helper methods in the
    parent class.  Furthermore, this design streamlines the ability of
    the manager to also operate as a worker when the baton assembler
    is being run as a single process.</p>

    <p>This class overrides some of the core methods defined in the
    Assembler (grandparent) class, which is the primary API for all
    assemblers in PEACE, thereby completing the class hierarchy that
    implements all the API methods.  Consequently, this class is
    instantiable.  However, this class should not and cannot be
    directly instantiated.  Instead, the AssemblerFactory class must
    be used to create a suitable instance of this class.  The
    AssemblerFactory is responsible for ensuring that the worker is
    instantiated only on processes whose MPI rank is greather than
    zero.</p>
*/
class BatonAssemblerWorker : public BatonAssembler {
    friend class AssemblerFactory;
public:
    /** The destructor.
        
        The destructor frees up all any dynamic memory allocated by
        this object for its operations.
    */
    virtual ~BatonAssemblerWorker();

    /** Primary method to perform cDNA assembly.
        
        This method is the primary interface method that is invoked
        (on each MPI process) to perform the complete assembly and
        output generation process.  This method is invoked once all
        the command line parameters are processed by the various
        subsystems constituting PEACE.  This method launches the
        various activities performed by the worker process.

		\note This method is invoked only after the initialize()
		method in the base class has successfully performed its
		operations.

        \return This method returns zero if the assembly process was
        successfully completed.  Otherwise this method returns a
        non-zero error code.
    */
    virtual int assemble();

protected:
	/** Helper method to send alignment information to the manager.

		This is a utility method that is used by the worker process to
		send alignment information to the manager.  This method is
		invoked once the worker has completed one round of checks
		against the subset of fragments assigned to it.  There are
		numerous calls to this method for each round of analysis
		initated by the manager.

		\note This method sends the alignment list to the manager as a
		flat list of characters.  This is a bit of a hack that we are
		exchanging structures as flat character arrays.  However, it
		works fine as long as we have a Single Program Multiple Data
		(SPMD) model and we are running on a homogeneous cluster

		\param[in] alignInfo The list of alignment information to be
		dispatched to the manager.  This list must always have a dummy
		entry at the end.  This last entry is ignored by the manager.
		However, having this dummy entry guarantees that every worker
		will always deliver an alignment list to the manager,
		streamlining the logic of handling distributed workers.
	*/
    void sendAlignmentInfo(const AlignmentInfoList& alignInfo) const;
    
private:
    /* The constructor for this class.
       
       The constructor has been made private so that this class cannot
       be directly instantiated.  However, notice that the
       AssemblerFactory is a friend of this class; therefore it can
       instantiate a BatonAssemblerWorker object.  Accordingly, the
       AssemblerFactory::create() method must be used to instantiate
       an object of this class.
    */
    BatonAssemblerWorker();
};

#endif
