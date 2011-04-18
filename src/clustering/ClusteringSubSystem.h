#ifndef CLUSTERING_SUB_SYSTEM_H
#define CLUSTERING_SUB_SYSTEM_H

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
#include "HeuristicChain.h"
#include "ArgParser.h"

// Forward declaration to keep compiler happy and fast
class ClusterMaker;
class ESTAnalyzer;

/** A major sub-system of PEACE that handles clustering of cDNA
    fragments.

	<p>As the name suggests, this class represents the filtering
	sub-system of PEACE.  It serves as the central point for
	performing various operations relating to clustering cDNA
	fragments.  The clustering is predominantly transcript-centric
	clustering where fragments from the same transcript (typically
	from the same gene) are placed in a single cluster.  This
	sub-system is run after the filtering sub-system.</p>

	<p>The clustering sub-system consists of a collection of
	ClusterMaker classes.  Each ClusterMaker performs a slightly
	different type of clustering. Each ClusterMaker uses an
	ESTAnalyzer object that is used to compare two cDNA fragments to
	determine a degree of similarity between them.  Not all
	ClusterMaker and ESTAnalyzer classes can be coupled together
	(nothing in the software prevents it but biologically it may not
	be meanigful to do so).</p>

    <p>Typically, the brunt of the clustering time is spent in the
	ESTAnalyzer to compare two cDNA fragments.  In order to improve
	performance each ESTAnalyzer is associated with a set of Heuristic
	classes organized in the form of a HeuristicChain.  Each Heuristic
	provides a different strategy to determine compatibility between
	two cDNA fragments, thereby helping the ESTAnalyzer to perform
	faster.  The FilterChain ties the heuristics together to yield
	more sophisticated heuristic mechanism to provide better tradeoff
	between performance and quality of the heuristics. Heuristics can
	be completely suppressed by specifying \c --heuristics \c null
	command-line argument.</p>
*/
class ClusteringSubSystem : public SubSystem {
public:
    /** The constructor.

        The constructor does not have very many tasks to perform other
        than to appropirately initialize member objects.  It is
        present to adhere to coding conventions and serve as a place
        holder for potential future enhancements.
    */
    ClusteringSubSystem();
    
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
    SubSystem::SubSystemKind getKind() const { return SubSystem::CLUSTERING; }

    /** Add the initial set of command line parameters for this
        sub-system.

        This method is invoked by the core system of PEACE just after
        the sub-system has been instantiated.  This method is expected
        to add the top-level command line arguments that it requires
        for its operation.  Accordingly, this method adds parameters
        for specifying the heuristics, EST analyzer, and cluster maker
        classes to be used for the current run of PEACE.
        
        \param[out] argParser The argument parser to which the basic
        command line arguments for this sub-system are to be added.
        This is the same command line parser that will be passed to
        the initialize() method.
    */
    void addCommandLineArguments(ArgParser& argParser);

    /** Method to initialize sub-system, create filters, and
        add them to a filter chain.

        This method is invoked after the initial round of command line
        processing has been completed by the core class.  As per the
        base class API description (see SubSystem::initialize() method
        documentation for general API contract for this method) this
        method performs the following tasks:

        <ol>

        <li>It creates the heuristics specified in the command line
        argument (ignoring the \c null heuristics name) and adds them
        to the HeuristicChain component in this subsystem.</li>

        <li>It then creates an valid ESTAnalyzer via the
        ESTAnalyzerFactory API method(s). It then arranges to have the
        EST analyzer's command-line parameters to be added to the
        argument parser (for second round of parameter
        processing)</li>

        <li>Next it creates the specified ClusterMaker via the
        ClusterMakerFactory API method(s).  It then arranges to have
        the Cluster Maker's command-line parameters to be added to the
        argument parser via the Component::addCommandLineArguments()
        method (for a second round of parameter processing)</li>

        </ol>

        \param[in,out] argParser The argument parser to which the
        additional command line arguments for this sub-system are to
        be added.  This is the same command line parser that was
        passed to the addCommandLineArguments() method.

        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur then this method must return an
        non-zero error code.
    */
    int initializeSubSystem(ArgParser& argParser);

    /** Initialize components in the sub-system.

        This method is invoked after the loadInputs() method has been
        invoked on all the sub-systems and the second round of command
        line parameter parsing has been completed.  As per the base
        class API description (see
        SubSystem::initializeSubComponents() method documentation for
        general API contract for this method) this method performs the
        following tasks:

        <ol>

        <li>It initializes the heurstic chain (causing the parameter
        set manager to be initialized).</li>

        <li>The ESTAnalyzer is initialized next.</li>

        <li>Finally, it initializes the cluster maker object.</li>

        </ol>
        
        \return This method returns NO_ERROR (zero) if no errors are
        encountered.  If errors occur then this method must return an
        non-zero error code.
    */
    virtual int initializeSubComponents();
    
    /** Invoked to permit the sub-system to perform its core tasks.

        This method is invoked once all the sub-systems have been
        successfully initialized.  This method performs the following
        tasks for the filtering sub-system:

        <ol>

        <li>First the heuristics in the chain are initialized via a
        call to the HeuristicChain::initialize() method.  That method
        call further initializes each Heuristic created and added to
        it in the initialize() method.</li>

        <li>If initialization of heuristics proceeds successfully,
        this method calls ESTAnalyzer::initialize that initializes the
        analyzer for this run of PEACE.</li>

        <li>Next, the ClusterMaker object for this run of PEACE is
        initialized.</li>

        <li>If errors occur during the above initialization steps
        (which includes validation of command-line parameters) then
        this method immediately exists with an non-zero value
        indicating error.</li>

        <li>Once all the components (and their sub-components) have
        been initialized (via the aforementioned steps) the core
        parallel clustering tasks are initiated via a call to
        ClusterMaker::makeClusters() method.</li>
        
        </ol>

        \return As per the API contract, this method currently returns
        zero if all the processing was successfully completed.
        On errors it returns a non-zero error code.
    */
    int run();

    /** Windup operations of this sub-system and its components.

        This method is invoked after all the sub-systems have
        successfully completed running.  In concordance with the API
        requirement for this method (see SubSystem::finalize() method
        documentation for general API contract for this method) this
        method directs each of the components to windup their
        operations.  Currently, this method merely invokes the
        corresponding finalization methods on each sub-system.
        
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
    void finalizeSubSystem(const bool success);
    
    /** The destructor.
        
        The destructor for this class.  Currently the destructor has
        no special tasks to perform (but is present to adhere to
        coding conventions).
    */
    virtual ~ClusteringSubSystem();
	
protected:
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has members
        whose value is set when the object is created.  These values
        cannot be changed during the lifetime of this object.
        
        \param[in] src The source object from where data is to be
        copied.  Currently this value is ignored.
        
        \return Reference to this.
    */
    ClusteringSubSystem& operator=(const ClusteringSubSystem& src);

    /** Component to handle the chain of heuristics specified by the
        user.

        This component is part of the ClusteringSubSystem. It deals
        with managing a list of Heuristic objects in the order
        specified by the user.
    */
    HeuristicChain heuristicChain;

    /** Reference to the ClusterMaker object being used for clustering.

        This instance variable is used to hold a pointer to the
        ClusterMaker object being used for clustering cDNA fragments
        for the current run of PEACE. This instance variable is
        initialized to NULL in the constructor.  The initialize()
        method uses the value in ClusteringSubSystem::clusterMakerName
        (filled-in by argument parsing) to instantiate a suitable
        derived ClusterMaker object via the ClusterMakerFactory. The
        finalize() method deletes this pointer and resets it back to
        NULL.
    */
    ClusterMaker *clusterMaker;

    /** Pointer to the ESTAnalyzer object being used for clustering.

        This instance variable is used to hold a pointer to the
        ESTAnalyzer object being used for comparing two cDNA fragments
        for the current run of PEACE. This instance variable is
        initialized to NULL in the constructor.  The initialize()
        method uses the value in ClusteringSubSystem::estAnalyzerName
        (filled-in by argument parsing) to instantiate a suitable
        derived ESTAnalyzer object via the ESTAnalyzerFactory. The
        finalize() method deletes this pointer and resets it back to
        NULL.
    */    
    ESTAnalyzer *estAnalyzer;
    
private:
    /** The command-line argument indicating the list and order in
        which heuristics are to be applied.

        This instance variable contains the command-line argument
        specified by the user for the \c --heuristics command-line
        parameter.  This list maintains the entries in the same order
        in which the user specified them. This list is used by the
        initialize() method to instantiate the corresponding Heuristic
        objects via the HeuristicFactory.
    */
    ArgParser::StringList heuristicNames;

    /** The command-line argument indicating the cluster maker to be used.

        This instance variable contains the command-line argument
        specified by the user for the \c --clusterMaker command-line
        parmaeter.  This string contains the name of the cluster maker
        to be used. This value is passed to the
        ClusterMakerFactory::create() method to instantiate a suitable
        cluster maker object.  The default value for this instance
        variable is \c "mst".
    */
    std::string clusterMakerName;

    /** The command-line argument indicating the EST analyzer to be
        used.

        This instance variable contains the command-line argument
        specified by the user for the \c --analyzer command-line
        parmaeter.  This string contains the name of the EST analyzer
        to be used. This value is passed to the
        ESTAnalyzerFactory::create() method to instantiate a suitable
        ESTAnalyzer object.  The default value for this instance
        variable is \c "twopassD2adapt".
    */    
    std::string estAnalyzerName;
};

#endif
