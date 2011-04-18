#ifndef CLUSTER_MAKER_H
#define CLUSTER_MAKER_H

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

#include "Component.h"
#include "ArgParser.h"

// Forward declaration to make compiler happy
class ESTAnalyzer;
class ESTList;
class MSTCluster;
class MST;

/** The base class of all cluster makers.

    This class must be the base class of all cluster makers in the
    system. This class provides some default functionality that can be
    readily used by each cluster maker.
*/
class ClusterMaker : public Component {
public:
   /** Add valid command line arguments for this analyzer.

        This method must be used to add all valid command line options
        that are supported by this analyzer.  Note that derived
        classes may override this method to add additional command
        line options that are applicable to it.  This method is
        invoked when the clustering sub-system is initialized.

        \note Derived EST analyzer classes may override this method to
        display help for their custom command line arguments.  When
        this method is overridden don't forget to call the
        corresponding base class implementation to add common options.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);
	
	/** A method to handle initialization tasks for the ClusterMaker.
		
		This method is called after the ESTs have been loaded and
		the clustering is just about to commence. This method
		calls the corresponding method in the base class to ensure
		base class initialization proceeds successfully.  It then
		sets up the ClusterMaker::estList convenience pointer in
		this class. Next, it inititalizes the ESTAnalyzer (which
		inturn initializes the HeuristicChain and the
		ParameterSetManager object in the HeuristicChain).
		
		\return This method returns \c true if initialization was
		successfully completed.  On errors it returns \c false.
	*/
    virtual bool initialize();
	
    /** Method to begin clustering.

        This method must be used to create clusters based on a given
        EST analysis method.  This method is a pure-virtual method.
        Therefore all cluster maker classes must override this method
        to perform all the necessary operations.

		\note This method must be invoked only after the initialize()
		method is invoked.
    */
    virtual int makeClusters() = 0;

    /** Obtain the EST analyzer set for this cluster maker.

        This method must be used to obtain the EST analyzer set for
        this cluster maker.

        \return The EST analyzer set for this cluster maker. If a
        valid EST analyzer has not been set then this method returns
        NULL.
    */
    inline ESTAnalyzer *getAnalyzer() const { return analyzer; }

    /** Add a dummy cluster to the cluster maker.

        This method can be used to add a dummy cluster to the cluster
        maker. The dummy clusters are added as direct descendants of
        the root cluster with the given name.

        \note This method is currently used by the Filter hierarchy to
        add ESTs that are logically filtered out and must not be part
        of the core clustering process.
        
        \param[in] name A human readable name to be set for this
        cluster.
        
        \return If a cluster was successfully added, then this method
        returns a unique integer that identifies the newly added
        cluster. This value must be used to add entries to this
        cluster via the ClusterMaker::addEST method.
    */
    virtual int addDummyCluster(const std::string name) = 0;

    /** Add a EST directly to a given cluster.

        This method can be used to add an EST directly to a
        cluster. This bypasses any traditional mechanism and directly
        adds the EST to the specified cluster.

        \note The EST is added with an invalid metric value. ESTs
        added to a cluster are not included in the standard clustering
        process. Adding an EST that has already been added to the
        same/another cluster results in undefined behaviors.

        \param[in] clusterID The unique ID of the cluster to which the
        EST is to be added. This value must have been obtained from an
        earlier (successful) call to the ClusterMaker::addDummyCluster
        method.
        
        \param[in] estIdx The EST to be added to the given
        cluster. Once the EST has been added to this cluster it will
        not be included in the clustering process performed by this
        cluster maker.
    */
    virtual void addEST(const int clusterID, const int estIdx) = 0;

	/** Obtain a pointer to the root cluster.

		This method can be used to obtain a reference to the set of
		clusters that has been built by this class.  This method
		returns a valid set of clusters only after the makeClusters()
		method has successfully completed.

		\return A pointer to the set of clusters built by this cluster
		maker.  The return pointer can be NULL.
	*/
	virtual const MSTCluster* getClusters() const = 0;

	/** Obtain a pointer to the Minimum Spanning Tree (MST) built by
		this class (if any).

		This method can be used to obtain an MST object built by this
		class (with help from one of the derived child classes).  This
		method returns a valid pointer only if the cluster maker
		builds an MST.  Otherwise this method returns NULL.

		\return The MST object built by this class.  The returned
		pointer can be NULL.
	*/
	virtual const MST* getMST() const = 0;
	
    /** The destructor.

        The destructor frees memory allocated for holding any data in
        this base class.
    */
    virtual ~ClusterMaker();

	
protected:
    /** The default constructor.

        The constructor has been made protected to ensure that this
        class is never directly instantiated.  Instead one of the
        derived cluster maker classes must be instantiated via the
        ClusterMakerFactor API methods.

        \param[in] name The human readable name for this cluster
        maker.  This name is used when generating errors, warnings,
        and other output messages for this object.

        \param[in,out] analyzer The EST analyzer to be used by this
        ClusterMaker for generating similarity metrics between two
        given ESTs.
    */
    ClusterMaker(const std::string& name, ESTAnalyzer *analyzer);

    /** The analyzer to be used for generating EST similarity metrics.

        This pointer is used to hold a pointer to the EST analyzer
        that must be used for generating similarity metrics between
        two given pairs of ESTs.  This pointer is initialized when the
        object is instantiated and is never changed during the
        lifetime of this object.
    */
    ESTAnalyzer* const analyzer;

    /** A shortcut reference to the shared list of cDNA fragments
        being analyzed.
		
        This list holds a pointer to the shared list of cDNA fragments
        currently being analyzed.  This pointer is initialized to NULL
        in the constructor.  A valid pointer is filled in when the
        initialize() method is invoked.  The pointer to the ESTList is
        obtained from the SubSystem::runtimeContext via the
        Component::subSystem member in the base class.
    */
    ESTList* estList;
	
private:
    /** A dummy operator=

        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this
        object.
        
        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.
        
        \return Reference to this.
    */
    ClusterMaker& operator=(const ClusterMaker& src);
};

#endif
