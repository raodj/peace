#ifndef SUB_SYSTEM_H
#define SUB_SYSTEM_H

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

// Forward declarations to keep compiler fast and happy
class ArgParser;
class RuntimeContext;

/** \file SubSystem.h

    This file defines the common base class for all primary
    sub-systems constituting PEACE.
*/

/** \def NO_ERROR 0

    \brief A named constant for zero.

    This named constant is used to remove magic numbers in the code
    and hopefully make the code a bit more readable.  Maybe this
    should become a strongly-typed static constant.
*/
#define NO_ERROR 0

/** A common, abstract base-class for a primary sub-system
    constituting PEACE.

    <p>PEACE has been designed as set of independent, but cooperative
    sub-systems that collectively provide its various features and
    functionality.  The objective is to achieve the desired level of
    "separation of concerns" so that each sub-system can independently
    evolve while seamlessly enabling interoperability between
    sub-systems.  For example, the design enables adding a new input
    or output file format without impacting the clustering or assembly
    algorithms.  Accordingly, the implementation strives to
    effectively use object-oriented design patterns and strategies to
    achieve the design objectives without sacrificing performance.</p>

    <p>This class serves a generic parent and interface class to every
    sub-system and establishes the most basic functionality that must
    be supported by every sub-system constituting PEACE.  Currently,
    PEACE is composed of the following five primary sub-systems: input
    sub-system, output sub-system, filtering sub-system, clustering
    sub-system, and assembly sub-system.  Note that each sub-system is
    implemented as class that extends this class and appropriately
    implements the required functionality.  In addition, each derived
    sub-system class provides mechanisms to create other
    sub-components and provides features specific to each sub-system.
    </p>

    <p>The general life cycle of the various sub-systems constituting
    PEACE is controlled by the primary PEACE class.  The life cycle
    proceeds in the following phases (after ll the sub-systems have
    been appropriately instantiated):

    <ol>

    <li><b>Phase 1</b>: The addCommandLineArguments() method is
    invoked on all the sub-systems (in some random order) to obtain
    the command line arguments supported by each one.</li>

	<li><b>Phase 2</b>: A first round of command-line argument parsing
    is performed by the top-level PEACE class.</li>

    <li><b>Phase 3</b>: Once the top-level command line arguments have
    been processed, the initializeSubSystem() method is invoked (in
    some random order) on each sub-system.  This method validates
    command line arguments and instantiates components.  Furthermore,
    additional command line parameters are gathered from the
    components.</li>

    <li><b>Phase 4</b>: A second round of command-line argument
    parsing is performed by the top-level PEACE class.  This round
    processes arguments associated with components.</li>

    <li><b>Phase 5</b>: The initializeComponents() method is invoked
    on all the sub-systems (in some random order).  This method can be
    optionally used to initialize components (however, input data set
    would not be available at this time). In addition, this method
    must establish valid pointers to shared components and data
    structures in the RuntimeContext.  <i>After this method call on
    all sub-systems, the shared RuntimeContext object will have all
    the necessary information filled-in with valid information based
    on the user-specified command-line arguments.</li>

	<li><b>Phase 6</b>: The loadInputs() method is invoked on all the
    sub-systems (in some random order).  This method must be used to
    load data files (after any pending initializations) into
    memory.</li>

    <li><b>Phase 7</b>: The initializeSubComponents() method is
    invoked on all the sub-systems (in some random order).  This
    method can be optionally used to further initialize components and
    any underlying sub-components (and primary input data sets would
    be available at this time).</li>
	
    <li><b>Phase 8</b>: The run() method is invoked on all the
    sub-systems in the following \e guaranteed order: input
    sub-system, output-subsystem, filtering sub-system, clustering
    sub-system, and assembly sub-system.  This method is where the
    core filtering, clustering, and assembly operations are
    performed.</li>

    <li><b>Phase 9</b>: The generateOutputs() method is invoked on all
    the sub-systems (in some random order).  This method is typically
    used to generate results from various operations performed by the
    different sub-systems.  Additional statistics about the run time
    performance and characteristics of the software sub-systems are
    also reported in this phase.</li>

    <li><b>Phase 10</b>: The finalizeComponents() method is invoked on
    all the subs-systems (in some random order).  This method is
    typically used to wind-down the operations of the various
    sub-components and components constituting the sub-systems.</li>

    <li><b>Phase 11</b>: The finalizeSubSystem() method is then
    invoked on all the sub-systems (in some undefined order).  This
    method can be used to wind-up the operations of the sub-system.
    </li>

    </ol>
    
    \note In the above life cycle, except for the run() method, there
    is no guarantees on the order in which each phase is invoked on
    the set of sub-systems.  That is, within a given phase, the order
    in which the appropriate life cycle API method is invoked is not
    guaranteed (and no attempt/assumptions must be made on such an
    order existing).  Any order that seems to exist is purely
    coincidental.

    </p>

    \note Since the life cycle is sufficiently elaborate (to
    accommodate various characteristics of the system and to have a
    good, symmetric API), this base class provide default
    implementations that perform absolutely nothing (other than return
    NO_ERROR where appropriate).  This is just a convenience so that
    derived classes don't have to implement all the phases of the
    overall system life cycle.  However, don't implement any
    functionality here.  This is meant to serve as a pure interface
    class.

*/
class SubSystem {
public:
    /** \brief Enumeration to identify the type of the sub-system.

        This enumeration list provides a convenient mechanism to
        identify the type of the sub-system represented by the derived
        class.
    */
    enum SubSystemKind {INVALID, CORE, INPUT, OUTPUT, FILTERING, CLUSTERING,
                        ASSEMBLY};

    /** Determine the type of sub-system represented by the derived
        class.

        This method can be used to determine the type of sub-system
        implemented by an object of this type. This method is
        implemented by the derived class to return an appropriate
        enumeration value.

        \return A valid enumeration constant indicating the type of
        the sub-system represented by this object.
    */
    virtual SubSystem::SubSystemKind getKind() const = 0;

    /** Add the initial set of command line parameters for this
        sub-system.

        This method is typically invoked by the core system of PEACE
        just after the sub-system has been instantiated.  This method
        is expected to add the top-level command line arguments that
        it requires for its operation. 

        \note The base class method implementation does absolutely
        nothing.
        
        \param[out] argParser The argument parser to which the basic
        command line arguments for this sub-system are to be added.
        This is the same command line parser that will be passed to
        the initialize() method.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

    /** Method to initialize sub-system, create components, and
        add additional arguments accepted by components.

        This method is invoked after the initial round of command line
        processing has been completed by the core class.  Each of the
        derived sub-systems are expected to perform the following
        tasks in this method:

        <ol>

        <li>First the components constituting the sub-system must be
        instantiated based on the command line argument values.</li>
        
        <li>Additional command-line arguments accepted by each
        component is then added to the supplied argument parser.
        These second-level arguments are parsed (by PEACE) after it
        collects the arguments from all sub-systems.<li>
        
        </ol>

        \note The method in this base class does absolutely nothing.
        
        \param[in,out] argParser The argument parser to which the
        additional command line arguments for this sub-system are to
        be added.  This is the same command line parser that was
        passed to the addCommandLineArguments() method.

        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
    */
    virtual int initializeSubSystem(ArgParser& argParser);

    /** Initialize components in the sub-system (if any).

        This method is invoked after the initializeSubSystem() method
        has been invoked on all the sub-systems and the second round
        of command line parameter parsing has been completed.

		\note At the time of this method call, the information
        regarding the set of cDNA fragments to be processed is not yet
        available.  If this information is needed for a component (or
        sub-component) initialization, then such tasks must be
        delegated to the initializeSubComponent() API call.
		
        Analogous to the initializeSubSystem() API call, this method
        is expected to perform the following tasks:

        <ol>

        <li>Initialize each component (and possibly sub-components) in
        the sub-system.  Initialization of components causes further
        validation of command line arguments specified by the
        user.</li>

        </ol>

        \note The default implementation in this base class does
        nothing and simply returns NO_ERROR.
        
        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur then this method must return an
        non-zero error code.
    */
    virtual int initializeComponents();
    
    /** Invoked to permit sub-systems to load any inputs they may
        need.
        
        This method is invoked after all the sub-systems and
        components have been successfully initialized.  This method
        must be used by the sub-system (or its components) to load
        necessary data into memory.  For example, in the input
        sub-system, this method would create appropriate file loaders
        and load data (or parts of data) into memory (taking possibly
        a few minutes).  Other sub-systems and components that may
        require an one-off input may read it when this method is
        invoked.

        \note The default implementation in this base class does
        nothing and simply returns NO_ERROR.
        
        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur during initialization then this
        method must return an non-zero error code.
    */
    virtual int loadInputs();

    /** Initialize sub-components in the sub-system (if any).

        This method is invoked after the loadInputs() method has been
        invoked on all the sub-systems.  At the time of this method
        call, the information regarding the set of cDNA fragments to
        be processed is ready and available.  This information is
        often needed for a component (or sub-component) initialization
        and such tasks can be performed by this method.  Analogous to
        the initializeSubSystem() API call, this method is expected to
        perform the following tasks:

        <ol>

        <li>Fully initialize each component and possibly
        sub-components in the sub-system.  This call cannot induce
        further changes to command line arguments.</li>
		
        </ol>
		
        \note The default implementation in this base class does
        nothing and simply returns NO_ERROR.
        
        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur then this method must return an
        non-zero error code.
    */
    virtual int initializeSubComponents();
	
    /** Invoked to permit the sub-system to perform its core tasks.
        
        This method is invoked once all the sub-systems and its
        components have been successfully initialized.  This method is
        expected to perform the following tasks:

        <ol>
        
        <li>It then arranges to have the components perform the core
        task of the sub-system.</li>
        
        </ol>

        For example, the clustering sub-system or assembly sub-system
        would perform the core clustering or assembly tasks (that may
        run for hours).  On the other hand, the input and output
        sub-systems may not have any operations to perform here.

        \note The default implementation in this base class does
        nothing and simply returns NO_ERROR.
        
        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur then this method must return an
        non-zero error code.
    */
    virtual int run();

    /** Generate outputs from this sub-system and its components.

        This method is invoked after all the sub-systems have
        successfully completed running.  This method is expected to
        generate any outputs from the sub-system.  This includes
        writing data to output files.  The output sub-system generates
        the data files during this method call.  Other sub-systems may
        generate outputs (such as statistics etc.) when this method is
        invoked.

        \note The default implementation in this base class does
        absolutely nothing.
        
        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    virtual void generateOutputs(const bool success);
    
    /** Windup operations of the components.

        This method is a symmetric-dual of the initializeComponents()
        method.  This method is invoked after all the sub-systems have
        completed running and outputs have been generated.  This
        method is expected to windup the operations of components.

        \note The default implementation in this base class does
        absolutely nothing.
        
        \param[in] success If this flag is \c true, then it indicates
        that the overall run of various sub-systems completed
        successfully.  Otherwise this flag is set to \c false to
        indicate at least one of the sub-systems reported an error.
    */
    virtual void finalizeComponents(const bool success);
    
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
        no special tasks to perform as this class is primarily just an
        interface class.  However, derived classes may perform some
        additional cleanup tasks here if needed.
    */
    virtual ~SubSystem();
    
    /** Set the runtime configuration that is shared by various
        sub-systems.

        This method is used by PEACE class to setup a pointer to the
        shared runtime configuration information.  The same runtime
        configuration information is shared between multiple
        sub-systems. This method simply saves the supplied pointer in
        the SubSystem::runtimeContext member for use by the derived
        class

        \param[in] context The runtime context/configuration
        information to be set for this sub-system.
    */
    void setContext(RuntimeContext* context);

    /** Obtain reference to the runtime context.

        This method must be used to obtain a reference to the shared
        runtime configuration information.

        \return The runtime context set for this sub-system.  If a
        valid runtime configuration has not been set, then this method
        returns NULL.
    */
    RuntimeContext* getContext() const { return runtimeContext; }

protected:
    /** The constructor.
        
        The only and default constructor for this class.  This class
        is not meant to be directly instantiated. Instead one of the
        derived classes must be appropriately instantiated.
    */
    SubSystem();

	/** The runtime context / configuration.

		This instance variable holds a pointer to the runtime context
		/ configuration information that is shared between multiple
		sub-systems.  This pointer is set by PEACE class when it
		creates various sub-systems.
	*/
    RuntimeContext* runtimeContext;
	
private:
    // Currently this class does not have any private members.
};

#endif
