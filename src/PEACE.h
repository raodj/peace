#ifndef PEACE_H
#define PEACE_H

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

#include "Utilities.h"
#include "RuntimeContext.h"
#include <vector>

// Forward declarations to avoid forcing users to include the various
// sub-systems.
class SubSystem;
class ArgParser;
class OutputSubSystem;
class InputSubSystem;
class FilteringSubSystem;
class ClusteringSubSystem;
class AssemblySubSystem;

/** The top-level class that controls and coordinates all the
    clustering and assembly activities.

    This is the top-level class constituting PEACE.  The primary
    method for triggering the various processing performed by this
    system is the run() method. The run() method (along with other
    helper methods in this class) perform the following core tasks:

    <ol>

	<li>The sub-systems constituting PEACE are instantiated when this
	class is created.</li>

	<li>It gathers the various top-level command line arguments that
	are accepted by various sub-systems by invoking the
	SubSystem::addCommandLineArguments() method.</li>
	
    <li>The command-line arguments supplied to PEACE (and passed-in
    via the global main() method) are processed for the initial,
    top-level options.  These options provide necessary information
    about various components of sub-systems to be used, such as:
    assembler to use, cluster generator to use, the set of filters to
    be applied, and the cDNA data files to be processed. Additional
    command-line arguments to the components of sub-systems are not
    processed and are deferred for processing by the components when
    the SubSystem::run() method is invoked.</li>

    <li>Assuming the top-level arguments parsed in the previous step
    are valid, the top-level components constituting the filtering,
    assembly, and clustering sub-systems are instantiated when the
    SubSystem::initialize() method is invoked on each sub-system. This
    method also adds command-line arguments that are specific to each
    component.  The sub-systems are initialized in the following
    order: <ol>

    <li>First, the output sub-system is initialized to create any
    debugging streams or logs (in one of many supported output
    formats) as indicated by the user.</li>

    <li>Next the file management sub-system is initialized to load any
    data files that the user may have specified.</li>
    
    <li>Next, the filtering sub-system is initialized.</li>
    
    <li>The cluster maker sub-system is initialized.</li>
    
    <li>Finally, the assembler sub-system is initialized.</li>
    
    </ol>
    
    </li>

    <li>Next, command-line arguments (based on parameters collected
    from each one of the sub-systems in the previous step during
    initialization) are processed and validated. The command-line
    arguments are processed in a 2-phased approach because the actual
    parameters vary depending on the type of sub-system(s) being
    utilized for a specific run of PEACE.</li>

    <li>If any errors occur in any of the aforementioned steps, then
    all valid command line arguments are displayed and no further
    processing is performed.</li>
    
    <li>Once all the sub-systems have been successfully initialized,
    they are all permitted to run in the same order as enumerated in
    the previous step.</li>
    
    <li>Finally, all the sub-systems are finalized in the following
    order:

	<ol>
        <li>The input sub-system is finalized. </li>

        <li>Next, the filtering sub-system is finalized.</li>

        <li>The cluster maker sub-system is finalized.</li>

        <li>The assembler sub-system is finalized.</li>

        <li>Lastly, the output sub-system is finalized closing all
        outputs.</li>
        		
        </ol>

	</li>
    
    </ol>
*/
class PEACE {
public:
    /** The constructor.

        The constructor does not have very many tasks to perform (as
        all the member objects appropirately initialize themselves).
        It establishes the shared runtime context information between
        the sub-systems.
    */
    PEACE();

    /** Method to initialize sub-systems, create sub-components, and
        add additional arguments accepted by sub-components.

        This method must be invoked before invoking the run() method.
        This method initializes all the sub-systems constituting
        PEACE. Initialization of sub-systems causes:

        <ul>

		<li> First the Message Passing Interface (MPI) library is
		initialized (assuming PEACE was compiled with MPI enabled)
		depending on the value of initMPI parameter.</li>

        <li> The command line arguments (if any) to be initially
        processed.</li>

        <li>Initialization of sub-systems (see
        SubSystem::initializeSubSystem() method for additional
        details) causes components to be instantiated.</li>

		<li>A second round of argument processing is performed and if
		any extraneous arguments remain, then this method displays
		usage and returns with an error.</li>

		<li>Next, each sub-system is directed to initialize its
		components via call to SubSystem::initializeComponents()
		method.</li>

		<li>Once the components have been initialized, this method
		calls subSystem::loadInputs() method to permit each sub-system
		and its components to load data required for their
		operation.</li>

		<li>Once the above method call(s) are successfully completed,
		this method invokes the SubSystem::initializeSubComponents()
		method on each sub-system.</li>
		
        </ul>

        \param[in,out] argc The number of command-line arguments that
        were passed-in to the process.  This instance variable
        determines the number of entries in the argv parameter.

        \param[in,out] argv The array of command-line arguments to be
        processed.

        \param[in] initMPI If this flag is \c true then this method
        initializes MPI.  Otherwise it is assumed that MPI
        initialization has been handled separately.  This enables
        multiple PEACE objects to be used in a single program run.
        
        \return This method returns zero on success.  If errors occur
        during initialization then this method must return an non-zero
        error code.
    */    
    int initialize(int& argc, char *argv[], bool initMPI = true);
    
    /** The core method that performs all the tasks.

        This method is the core method that is called the global
        main() function to perform the tasks associated with the
        current invocation of PEACE.

        \note If the PEACE::initialize() method has not be called
        before calling run() then this method calls
        PEACE::initialize() first.

        \param[in,out] argc The number of command-line arguments that
        were passed-in to the process.  This instance variable
        determines the number of entries in the argv parameter.

        \param[in,out] argv The array of command-line arguments to be
        processed.

        \param[in] initMPI If this flag is \c true then this method
        initializes MPI.  Otherwise it is assumed that MPI
        initialization has been handled separately.  This enables
        multiple PEACE objects to be used in a single program run.
        
        \return This method returns an exit code indicating the
        overall success status of invoking this method. On success,
        this method returns zero.  Otherwise it returns a non-zero
        value indicating error.
    */
    int run(int& argc, char *argv[], bool initMPI = true);

	/** Performs all clean-up operations as needed.

		This method must be called after the run() method has been
		invoked.  This method performs clean-up operations on the
		sub-systems and their components.  This method optionally
		finalizes MPI depending on the value of finalizeMPI parameter.

        \param[in] finalizeMPI If this flag is \c true then this
        method finalizes MPI.  Otherwise it is assumed that MPI
        finalization would be handled separately.  This enables
        multiple PEACE objects to be used in a single program run.

		\note This method must be called immaterial of the exit status
		from the run() method.  The destructor also calls this method
		if needed.
	*/
	void finalize(bool finalizeMPI = true);
	
    /** The destructor.
        
		The destructor checks and calls the finalize() method (if it
        was not called as required by the API).
    */
    virtual ~PEACE();

	/** Obtain a mutable pointer to the shared runtime context.

		This method can be used to obtain a pointer to the shared
		runtime context object associated with this instance of PEACE.

		\return This method returns a valid pointer to the current
		runtime context.  This pointer is never NULL.
	*/
	inline RuntimeContext* getContext() { return &runtimeContext; }

	/** Obtain an immutable pointer to the shared runtime context.

		This method can be used to obtain a constant pointer to the
		shared runtime context object associated with this instance of
		PEACE.

		\return This method returns a valid constant pointer to the
		current runtime context.  This pointer is never NULL.
	*/
	inline const RuntimeContext* getContext() const { return &runtimeContext; }
	
    /** Helper method to load a given file with a given optional file
        format.

        This method is a helper method that can be used to load cDNA
        fragments from a given data file. This method loads cDNA
        fragments from the given file into the list of ESTs maintained
        in the runtime context (RuntimeContext::estList).

        \note This method uses the InputFileFactory in the
        InputSubSystem to load the data file.
        
        \param[in] fileName The file (with optional relative or
        absolute path) to be loaded by this method.

        \param[in] format The format of the data file to be loaded.
        Valid strings are \c "fasta" or \c "sff".  By default the file
        type is auto detected and loaded.

        \param[in] startIndex The starting index position of the EST
	from where the full data is to be retained. For other entries,
	the header information and nucleotide sequences are loaded
	on-demand.  This is done to reduce peak memory footprint
	(thereby enabling the processing of large data files). This
	value must be less than endIndex.  This value is with
	reference to the full list of ESTs (not just relative to one
	file).
        
        \param[in] endIndex The ending index position of the EST upto
	(and <b>not</b> including) which the full data is to be
	retained. For other entries, the header information and
	nucleotide sequences are loaded on-demand.  This value must be
	greater than startIndex. This value is with-reference-to the
	full list of cDNA fragments.  By default the startIndex and
	endIndex are set such that all cDNA fragments are loaded
	on-demand.
        
        \return If all the files is successfully loaded, then this
        method returns \c true.  If an error occurs when processing
        one the files, this method returns with \c false.
    */    
    bool loadFile(const std::string& fileName, const std::string& format = "",
                  const long startIndex = MAX_READS,
                  const long endIndex   = MAX_READS);

protected:
    /** The output sub-system.

        This member object serves as the output sub-system for a given
        instance of PEACE.  It serves as the central point for
        performing various output operations from the other
        sub-systems constituting PEACE.  It is used for generating
        output files in different file formats, debugging outputs, log
        messages, and error messages.

		\note The object is created in the constructor and deleted in
		the destructor.
    */
    OutputSubSystem* oss;

	/** The input sub-system.
		
        This member object serves as the input sub-system for a given
        instance of PEACE.  It serves as the central point for various
        input operations.  It loads data (stored in different formats)
        used by the other sub-systems constituting PEACE.

		\note The object is created in the constructor and deleted in
		the destructor.		
	*/
    InputSubSystem* iss;

    /** The filtering sub-system.
		
        This member object serves as the filtering sub-system for a
        given instance of PEACE.  It serves as the central point for
        various filtering operations.  It essentially filters out cDNA
        fragments that could potentially degrade the quality of
        clustering/assembly.

		\note The object is created in the constructor and deleted in
		the destructor.
    */
    FilteringSubSystem* fss;
	
	/** The clustering sub-system.
        
        This member object serves as the clustering sub-system for a
        given instance of PEACE.  It serves as the central point for
        various clustering operations.  It essentially places a set of
        related (arising from same transcript or gene) cDNA fragments
        into a single cluster.

		\note The object is created in the constructor and deleted in
		the destructor.
    */
    ClusteringSubSystem* css;

	/** The assembly sub-system.
        
        This member object serves as the assembly sub-system for a
        given instance of PEACE.  It serves as the central point for
        various cDNA assembly operations.  It essentially places a set
        of related (arising from same transcript or gene) cDNA
        fragments into a single contig.

		\note The object is created in the constructor and deleted in
		the destructor.
    */
    AssemblySubSystem* ayss;
	
private:
	/** The shared runtime context information.

		This member maintains a simple hash map that is used to share
		runtime configuration information between various sub-systems
		constituting PEACE.  The constructor shares this runtime
		configuration object between multiple systems via
		SubSystem::setContext() method.
	*/
	RuntimeContext runtimeContext;

    /** A temporary argument parser shared between initialize() and
        run() methods.

        This is a temporary argument parser that is created in the
        initialize() method and deleted in the finalize()
        method. Although its use is rather short, it is used to ensure
        that the various API methods are called to ensure consistent
        startup and shutdown.  Consequently, it has been defined as a
        pointer and not as an object.
    */
    ArgParser *argParser;

	/** Instance variable to track the last error code.

		This instance variable is used to track the last error code
		that was reported by any method in this class.  This value is
		initialized to zero.  When an error occurs, it is set to a
		non-zero value (and never changed).
	*/
	int errorCode;

    /** A convenience vector to which sub-systems are added and
        processed as a batch.

        The constructor instantiates and adds pointers to various
        sub-systems to this vector.  The various operations that are
        to be performed on the sub-systems are achieved by iterating
        over this list and invoking appropriate methods on them.
    */
    std::vector<SubSystem*> subSystemList;
};

#endif
