#ifndef INPUT_SUB_SYSTEM_H
#define INPUT_SUB_SYSTEM_H

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
#include "InputFileFactory.h"
#include "OnDemandESTList.h"

/** A major sub-system of PEACE that handles loading of different data
    from various data sources.

	<p>As the name suggests, this class represents the input
	sub-system of PEACE.  It serves as the central point for
	performing various input operations to load data required or used
	by the other sub-systems constituting PEACE.  It is used for
	reading cDNA fragments (to be processed) from different file
	formats, reading intermediate results (for reusing previous
	computations), and configuration data.</p>

	<p>Specifically, this sub-system consists of the following
	components:

	<ol>

	<li>A InputFileFactory component that is capable of creating and
	loading cDNA fragments from different file formats.  Each type of
	file format is read using sub-components derived from the
	InputFile class.</li>

	<li>A MSTReader component that can load Minimum Spanning Tree
	(MST) data that was previously written by PEACE.</li>

	</ol>
		
	</p>
*/
class InputSubSystem : public SubSystem {
public:
    /** The constructor.

        The constructor does not have very many tasks to perform other
        than to appropirately initialize member objects.  It is
        present to adhere to coding conventions and serve as a place
        holder for potential future enhancements.
    */
    InputSubSystem();
    
    /** Determine the type of sub-system represented by the derived
        class.

        This method can be used to determine the type of sub-system
        implemented by an object of this type. This method overrides
        the pure-virtual declaration in the base-class as mandated by
        the base class API.

        \return A valid enumeration constant indicating the type of
        the sub-system represented by this object.  This method always
        returns SubSystem::INPUT.
    */
    SubSystem::SubSystemKind getKind() const { return SubSystem::INPUT; }

    /** Add the initial set of command line parameters for this
        sub-system.

        This method is invoked by the core system of PEACE just after
        the sub-system has been instantiated.  This method is expected
        to add the top-level command line arguments that it requires
        for its operation.  Accordingly, this method adds parameters
        for various input components, such as: the cDNA file to be
        loaded, options to handle masked bases and \c N nucleotide
        entries etc. Most of these options are actually added by the
        InputFileFactory component and the MSTReader
        component.
        
        \param[out] argParser The argument parser to which the basic
        command line arguments for this sub-system are to be added.
        This is the same command line parser that will be passed to
        the initialize() method.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Method to initialize sub-system, create components, and
        add additional arguments accepted by components.

        This method is invoked after the initial round of command line
        processing has been completed by the core class.  As per the
        base class API description (see
        SubSystem::initializeSubSystem() method documentation for
        general API contract for this method) this method performs the
        following tasks:

        <ol>

        <li>First it creates a suitable instance of the ESTList object
        (that will hold the set of cDNA objects to be processed) and
        sets up the reference in the SubSystem::runtimeContext object
        (which is shared by each SubSystem constituting PEACE).</li>
        
        <li>Then it creates the InputFileFactory component.</li>
        initialize and load the cDNA data from the data files
        specified by the user (if any).</li>
        
        <li>Next, it checks and creates the MSTInput component which
        loads MST data from the file specified by the user (if
        any).</li>

        <li>Next, the ClusterInput component is instantiated.  This
        component is used to load clustering results from the file
        specified by the user (if any).</li>

        </ol>

        \param[in,out] argParser The argument parser to which the
        additional command line arguments for this sub-system are to
        be added.  This is the same command line parser that was
        passed to the initialize() method.

        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
    */
    int initializeSubSystem(ArgParser& argParser);

    /** Invoked to permit sub-systems to load any inputs they may
        need.
        
        This method is invoked after the initialiseSubSystem() method
        has been called. As per the base class API description (see
        SubSystem::loadinputs() method documentation for general API
        contract for this method) this method performs the following
        tasks:

        <ol>

        <li>First, the list of cDNA fragments are processed and
        indexed (for full loading on demand) into the shared
        RuntimeContext::estList.</li>
        
        <li>It directs the MSTInput component to load MST data from
        the file specified by the user (if any).</li>

        <li>Next, the ClusterInput component is initialized to load
        clustering results from the file specified by the user (if
        any).</li>

        </ol>
        
        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
	*/
	int loadInputs();

	/** Obtain a reference to the input file file factory.

		\return A reference to the input file factory component
		associated with this sub-system.
	 */
	InputFileFactory& getInputFileFactory() { return inputFileFactory; }
	
    /** Windup operations of the components.

        This method is invoked after all the sub-systems have
        successfully completed running and outputs have been
        generated.  In concordance with the API requirement for this
        method (see SubSystem::finalizeComponents() method
        documentation for general API contract for this method) this
        method directs each of the components to windup their
        operations and close all output streams/files.

        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    void finalizeComponents(const bool success);
    
    /** The destructor.
        
        The destructor for this class.  Currently the destructor has
        no special tasks to perform (but is present to adhere to
        coding conventions).
    */
    virtual ~InputSubSystem();
	
protected:
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has members
        whose value is set when the object is created.  These values
        cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    InputSubSystem& operator=(const InputSubSystem& src);
    
private:
    /** Component to handle loading of cDNA entries from various data files.

        This component is part of the InputSubSystem. It deals with
        loading cDNA entries from data files with different file
        formats (such as: FASTA, SFF, SAM, BAM etc.).
    */
    InputFileFactory inputFileFactory;

    /** The shared EST list.

        This instance variable contains the shared EST list that is
        used by various sub-systems. A pointer to the shared EST list
        is stored in the RuntimeContext::estList and 
     */
    OnDemandESTList sharedESTList;
};

#endif
