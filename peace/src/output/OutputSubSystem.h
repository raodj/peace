#ifndef OUTPUT_SUB_SYSTEM_H
#define OUTPUT_SUB_SYSTEM_H

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
#include "StandardOutput.h"
#include "MSTWriter.h"
#include "ClusterWriter.h"

/** A major sub-system of PEACE that handles generation of various
    outputs from components constituting PEACE.

	<p>As the name suggests, this class represents the output
	sub-system of PEACE.  It serves as the central point for
	performing various output operations from the other sub-systems
	constituting PEACE.  It is used for generating output files in
	different file formats, debugging outputs, log messages, and error
	messages.  Some of these output file formats can be customized
	using suitable command line parameters.</p>

	<p>Specifically, this sub-system consists of the following
	components:

	<ol>

	<li>A general purpose output stream that is the same as std::cout.
	However, the output stream could be redirected to a file through
	the command-line argument (\c --stdout).</li>

	<li>A logger stream that is tied to std::clog and can be used to
	generate logging message.  However, the output stream could be
	redirected to a file through the command-line argument (\c
	--stdclog).  In addition, the sub-system modifies the logger
	stream such that it logs only a subset of messages based on the
	current log level set by the user (\c --logLevel=info).</li>

	<li>A specialized component to write clustering solution to a
	given destination.  The cluster output component provides several
	different sub-components that can generate the final output in a
	variety of different formats that are suitable for other systems
	(such as: PEACE-GUI).</li>

	<li>A specialized component to write the Minimum Spanning Tree
	(MST) data structure to a given destination.</li>

	</ol>
		
	</p>
*/
class OutputSubSystem : public SubSystem {
public:
    /** The constructor.

        The constructor does not have very many tasks to perform other
        than to appropirately initialize member objects.  It is
        present to adhere to coding conventions and serve as a place
        holder for potential future enhancements.
    */
    OutputSubSystem();
    
    /** Determine the type of sub-system represented by the derived
        class.

        This method can be used to determine the type of sub-system
        implemented by an object of this type. This method overrides
        the pure-virtual declaration in the base-class as mandated by
        the base class API.

        \return A valid enumeration constant indicating the type of
        the sub-system represented by this object.  This method always
        returns SubSystem::OUTPUT.
    */
    SubSystem::SubSystemKind getKind() const { return SubSystem::OUTPUT; }

    /** Add the initial set of command line parameters for this
        sub-system.

        This method is invoked by the core system of PEACE just after
        the sub-system has been instantiated.  This method is expected
        to add the top-level command line arguments that it requires
        for its operation.  Accordingly, this method adds parameters
        for various output components (such as: standard output,
        standard error, standard log, MSTOutput, ClusterOutput)
        supported by this sub-system.
        
        \param[out] argParser The argument parser to which the basic
        command line arguments for this sub-system are to be added.
        This is the same command line parser that will be passed to
        the initialize() method.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Method to initialize components and add additional arguments
        accepted by components.

        This method is invoked after the initial round of command line
        processing has been completed by the core class.  As per the
        base class API description (see
        SubSystem::initializeComponents() method documentation for
        general API contract for this method) this method performs the
        following tasks:

        <ol>

        <li>First it permits the StandardOutput component to suitably
        handle any additional redirection or configuration of the
        standard output streams (namely: std::cout and
        std::clog).</li>
        
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

    /** Method to initialize components and add additional arguments
        accepted by components.

        This method is invoked after the second round of command line
        processing has been completed by the core class.  As per the
        base class API description (see
        SubSystem::initializeComponents() method documentation for
        general API contract for this method) this method performs the
        following tasks:

        <ol>

        <li>It permits the MSTOutput component is initialized to
        validate and create an appropriate output sub-component that
        can serialize a Minimum Spanning Tree (MST) data structure.</li>

        <li>The ClusterOutput component is initialized next to create
        an appropriate output sub-component that serializes the
        clustering output in an appropriate format.</li>

        <li>Lastly, the AssemblyOutput component is initialized to
        create an appropriate sub-component to generate assembly
        output in different file formats.</li>
        
        </ol>

        \param[in,out] argParser The argument parser to which the
        additional command line arguments for this sub-system are to
        be added.  This is the same command line parser that was
        passed to the initialize() method.

        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
    */
    int initializeComponents(ArgParser& argParser);

    /** Generate outputs from this sub-system and its components.

        This method is invoked after all the sub-systems have
        successfully completed running.  This method is expected to
        generate any outputs from the sub-system.  This includes
        writing data to output files.  The output sub-system generates
        the data files during this method call.
        
        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    virtual void generateOutputs(const bool success);
    
    /** Windup operations of this sub-system.

        This method is invoked after all the sub-systems have
        successfully completed running.  In concordance with the API
        requirement for this method (see
        SubSystem::finalizeSubSystem() method documentation for
        general API contract for this method) this method directs each
        of the components to windup their operations and close all
        output streams/files.

        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    void finalizeSubSystem(const bool success);
    
    /** The destructor.
        
        The destructor for this class.  Currently the destructor has
        no special tasks to perform (but is present to adhere to
        coding conventions).
    */
    virtual ~OutputSubSystem();

	/** Obtain the MST writer to be used for serializing a Minimum
		Spanning Tree.

		This method must be used to obtain a reference to the
		component to be used for serializing a Minimum Spanning Tree
		(MST) data to a file.

		\return The MST writer component associated with this
		sub-system.
	*/
	MSTWriter& getMSTWriter() { return mstWriter; }
	
protected:
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has members
        whose value is set when the object is created.  These values
        cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    OutputSubSystem& operator=(const OutputSubSystem& src);
    
private:
    /** Component to redirect standard output streams.

        This component is part of the OutputSubSystem. It deals with
        redirecting standard output (std::cout) and standard log
        (std::clog) to files.
    */
    StandardOutput stdOutput;

	/** Component to write MST data.

		This component is part of the OutputSubSystem. It deals with
		serializing the MST data structure to a given file.
	*/
	MSTWriter mstWriter;

	/** Component to write clustering data.

		This component is part of the OutputSubSystem. It deals with
		serializing the clustering data structure (represented by the
		MSTCluster class) to a given file.
	*/
    ClusterWriter clusterWriter;
    
    /** The command-line argument indicating the output FASTA file to
        which ESTs that where filtered out must be written.

        This instance variable contains the name of the file
        (specified by the user as a command-line argument to \c
        --filter-fail) to which cDNA fragments that were filtered out
        must be written.  The entries are written in a FASTA format.
    */
    std::string filterFailFile;

    /** The command-line argument indicating the output FASTA file to
        which ESTs that where not filtered-out must be written.

        This instance variable contains the name of the file
        (specified by the user as a command-line argument to \c
        --filter-pass) to which cDNA fragments that were <b>not</b>
        filtered out must be written.  The entries are written in a
        FASTA format.
    */
    std::string filterPassFile;	
};

#endif
