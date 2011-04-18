#ifndef ASSEMBLY_SUB_SYSTEM_H
#define ASSEMBLY_SUB_SYSTEM_H

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

#include "SubSystem.h"
#include <string>

// Forward declaration to keep compiler happy and fast
class ClusterMaker;
class ESTAnalyzer;
class ArgParser;

/** A major sub-system of PEACE that handles assembly of cDNA
    fragments into contigs.

	<p>As the name suggests, this class represents the assembly
	sub-system of PEACE.  It serves as the central point for
	performing various operations relating to assembly of cDNA
	fragments.  Assembly process can proceed using one of the
	following two approaches:

    <ul>

    <li>In one approach, an explicit phase of clustering is performed
	to cluster sufficiently similar cDNA fragments into a single
	cluster. In the next phase the fragments in each cluster are
	analyzed to form one (or more) contigs.</li>

    <li>Alternatively, the process of clustering and contig formation
    is blended together and occur simultaneously.  Consequently,
    unlike the former approach, there isn't a very clear dichotomy
    between two phases</li>
    
    </ul>
    
	<p>The assembly sub-system consists of a collection of Assembler
	classes.  Each Assembler performs a different type of
	assembly. Each Assembler may optionally use a suitable ESTAnalyzer
	object that is used to compare two cDNA fragments to determine a
	degree of similarity between them.  Not all Assembler and
	ESTAnalyzer classes can be coupled together (nothing in the
	software prevents it but algorithmically it may not be meanigful
	to do so).</p>
*/
class AssemblySubSystem : public SubSystem {
public:
    /** The constructor.

        The constructor does not have very many tasks to perform other
        than to appropirately initialize member objects.  It is
        present to adhere to coding conventions and serve as a place
        holder for potential future enhancements.
    */
    AssemblySubSystem();
    
    /** Determine the type of sub-system represented by the derived
        class.

        This method can be used to determine the type of sub-system
        implemented by an object of this type. This method overrides
        the pure-virtual declaration in the base-class as mandated by
        the base class API.

        \return A valid enumeration constant indicating the type of
        the sub-system represented by this object.  This method always
        returns SubSystem::CLUSTERING.
    */
    SubSystem::SubSystemKind getKind() const { return SubSystem::ASSEMBLY; }

    /** Add the initial set of command line parameters for this
        sub-system.

        This method is invoked by the core system of PEACE just after
        the sub-system has been instantiated.  This method is expected
        to add the top-level command line arguments that it requires
        for its operation.  Accordingly, this method adds parameters
        for specifying the assembler (and optionally an EST
        analyzer/cluster maker) to be used for the current run of
        PEACE.
        
        \param[out] argParser The argument parser to which the basic
        command line arguments for this sub-system are to be added.
        This is the same command line parser that will be passed to
        the initialize() method.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Method to initialize sub-system and create an assembler.

        This method is invoked after the initial round of command line
        processing has been completed by the core class.  As per the
        base class API description (see
        SubSystem::initializeSubSystem() method documentation for
        general API contract for this method) this method performs the
        following tasks:

        <ol>

        <li>It then creates an valid assembler object via the
        AssemblerFactory API method(s). It then arranges to have the
        assembler's command-line parameters to be added to the
        argument parser (for second round of parameter
        processing).</li>

        <li>It sets up the shared pointer to the assembler in the
        RuntimeContext::assembler instance variable.</li>
        
        </ol>

        \param[in,out] argParser The argument parser to which the
        additional command line arguments for this sub-system are to
        be added.  This is the same command line parser that was
        passed to the addCommandLineArguments() method.

        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
    */
    int initializeSubSystem(ArgParser& argParser);

    /** Method to initialize components and sub-components.

        This method is invoked after the second round of command line
        processing has been completed by the core class.  As per the
        base class API description (see
        SubSystem::initializeSubComponents() method documentation for
        general API contract for this method) this method performs the
        following tasks:

        <ol>

        <li>It initializes the assembler object which validates its
        command-line parameters etc.</li>

        </ol>

        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
    */    
    int initializeSubComponents();

    /** Invoked to permit the sub-system to perform its core tasks.

        This method is invoked once all the sub-systems have been
        successfully initialized.  This method performs the following
        tasks for the filtering sub-system:

        <ol>

        <li>If a valid assembler object is available, it triggers the
        assembly operation via call to Assembler::assemble()
        method.</li>
        
        </ol>
        
        \return As per the API contract, this method currently returns
        zero if all the processing was successfully completed.
        On errors it returns a non-zero error code.
    */
    int run();

    /** Windup operations of the components in the sub-system.

        This method is invoked after all the sub-systems have
        successfully completed running.  In concordance with the API
        requirement for this method (see
        SubSystem::finalizeComponents() method documentation for
        general API contract for this method) this method directs each
        of the components to windup their operations.  Currently, this
        method merely invokes the corresponding finalization methods
        on each component.
        
        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    void finalizeComponents(const bool success);

    /** Windup operations of this sub-system as a whole.

        This method is a symmetric-dual of the initializeSubSystem()
        method.  This method is invoked after the finalizeComponents()
        method has been invoked on all the sub-systems.  This method
        is expected to windup the operations of the sub-system by
        deleting components (objects) and freeing up dynamic memory
        associated with the sub-system.

        \note The default implementation in this base class does
        absolutely nothing.
        
        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    virtual void finalizeSubSystem(const bool success);
    
    /** The destructor.
        
        The destructor for this class.  Currently the destructor has
        no special tasks to perform (but is present to adhere to
        coding conventions).
    */
    virtual ~AssemblySubSystem();
	
protected:
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has members
        whose value is set when the object is created.  These values
        cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    AssemblySubSystem& operator=(const AssemblySubSystem& src);

private:
    /** The command-line argument indicating the assembler to be used.

        This instance variable contains the command-line argument
        specified by the user for the \c --assembler command-line
        parmaeter.  This string contains the name of the assembler
        This value is passed to the AssemblerFactory::create() method
        to instantiate a suitable assembler object.  This variable
        does not have a valid default value.
    */    
    std::string assemblerName;
};

#endif
