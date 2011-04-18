#ifndef MST_CLUSTER_H
#define MST_CLUSTER_H

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

#include "MST.h"
#include "ESTAnalyzer.h"

class MSTCluster;
class ESTList;

/** \def ClusterList

    \brief Typedef for std::vector<MSTCluster>.

    This typedef is a convenience definition to handle a vector of
    Cluster objects associated with a given MSTCluster.
*/
typedef std::vector<MSTCluster*> ClusterList;

/** Class to encapsulate cluster information and aid in building
    clusters using Minimum Spanning Tree (MST) data.

    This class was introduced to encapsulate and manage the data
    associated with various clusters built using MST information.
*/
class MSTCluster {
    // Insertion operator to make dumping MST for debugging easier.
    friend std::ostream& operator<<(std::ostream&, const MSTCluster&);  
public:
    MSTCluster(MSTCluster* owner = NULL, const std::string& name = "");
    ~MSTCluster();

    /** Method to make clusters.

        This method takes a supplied list of nodes from an MST,
        an EST analyzer, and a threshold score, and uses them to
        break the MST into clusters based on the threshold.

        \param nodeList  The list of nodes in the MST.
        \param analyzer  The ESTAnalyzer used for score computation.
        \param threshold The threshold value for separating clusters.
    */
    void makeClusters(const NodeList& nodeList, const ESTAnalyzer* analyzer,
                      const float threshold);

    /** Method to automatically calculate a threshold.

        Currently this method is called only when the threshold supplied
        to MSTCluster::makeClusters() is equal to -1.  This threshold
        corresponds to the CLU analyzer and tells the calculateThreshold
        method to calculate the threshold based on the mean and variance
        of scoring metrics among the list of MST nodes.

        \param nodeList The list of nodes in the MST

        \return the calculated threshold as a floating-point number
    */
    float calculateThreshold(const NodeList& nodeList) const;

    /** Add a MST Node to this cluster.

        \param[in] node The node ot be added to the cluster.
    */
    void add(const MSTNode& node);

    /** Add a child cluster to this cluster.

        This method must be used to add a child cluster this
        cluster. Note that the parent cluster must be appropriately
        set in the child to refer to this MSTCluster.

        \note This cluster node takes ownership of the child node.

        \param[in] child The child node to be added to the this
        cluster.
    */
    void add(MSTCluster *child);
    
    /** Method to directly add a given EST to a dummy clusterID.

        This method is typically used to add a EST that has been
        filtered out by a filter to a specific cluster.  The cluster
        maker typically creates a dummy MSTNode with the necessary
        information to be added.
        
        \param clusterID The ID/index of the cluster to which the
        MST node is to be added.
        
        \param[in] node The MSTNode entry to be added to this cluster.
    */
    void add(const int clusterID, const MSTNode& node);

    /** Print the cluster tree.

        This method prints a tree display of the clusters.
    */
    void printClusterTree(const ESTList& estList,
                          std::ostream& os = std::cout,
                          const std::string& prefix = "") const;

    /** Print the cluster tree in a different format.

        As above, but prints the cluster tree in a format recognizable
        by the PEACE GUI and suitable for GUI processing.
    */
    void guiPrintClusterTree(std::ostream& os = std::cout) const;

    /** makeMergedClusters

        Currently not used in the code.
    */
    void makeMergedClusters(const int size, int* parent, bool* root);

    /** Gets the ID of this cluster.

        \return the cluster ID, an integer
    */
    inline int getClusterID() const { return clusterID; }

	/** Determine if this MSTCluster node has some members.

		This method can be used to determine if this cluster node has
		any sub-members in it. It considers both sub-clusters and cDNA
		fragments contained by it.

		\return This method returns \c false if the node has some
		members.  Otherwise this method returns \c true indicating it
		is empty.
	*/
	inline bool isEmpty() const {
	  return members.empty() && clusterList.empty();
	}
		
protected:
    /** A helper method for the guiPrintClusterTree method.

        Works by printing the data for the current cluster
        and then calling guiPrintTree on the sub-clusters (if any).
    */
    void guiPrintTree(std::ostream& os) const;

	
    /** Returns the FASTA header corresponding to a given EST.
        
		This method is a helper method that is used to format the
		information associated with a given EST for display purposes.
		The necessary information about the EST is determined via the
		EST::getInfo() method.

        \param[in] estIdx The index of the EST in the estList whose
        information is to be formatted and returned by this method.

        \param[in] estList The list of ESTs from where the information
        is to be extracted.
        
		\return The FASTA header associated with the given EST.  If no
		information is available this method returns an empty string.
    */
	std::string getESTInfo(const int estIdx, const ESTList& estList) const;
    
private:
    /** List of sub-clusters for this cluster.

        This instance variable is used to maintain the list of
        sub-clusters for this cluster. Entries are added to this list
        whenever sub-clusters are created in the
        MSTCluster::makeClusters() method.
    */
    ClusterList clusterList;
    
    const MSTCluster* parent;
    NodeList members;
    const int clusterID;

    /** A name set to identify filtered clusters.

        The name is set when dummy clusters are created to add ESTs
        that were filtered out based on a specific condition.  The
        named clusters are typically created by filters.  By default
        clusters don't have a name.  These indicate regular clusters.
        
        \see Filter.
    */
    const std::string name;
	
    /** Instance variable to track the next available cluster ID.
        
        This instance variable is used to generate unique cluster ID
        values for each newly created cluster.  It is intialized to
        zero. Each time a cluster is instantiated, the constructor
        uses this value to set the clusterID and then increments this
        value.
    */
    static int clusterIDSequence;

    /** Global list to maintain reference to all the MSTClusters created.

        This list is used to maintain a pointer to all the MST
        clusters ever created. This list is used to look up clusters
        given the index of the cluster.
		
        \see MSTCluster::getCluster() method
    */
    static ClusterList globalClusterList;
    
    /** A dummy operator=
        
        The operator=() is supressed for this class as it has constant
        members whose value is set when the object is created.  These
        values cannot be changed during the lifetime of this object.

        \param[in] src The source object from where data is to be copied.
        Currently this value is ignored.
        
        \return Reference to this.
    */
    MSTCluster& operator=(const MSTCluster& src);
};

/** \fn std::ostream& operator<<(std::ostream&, const MSTCluster&)

    Insertion operator to stream MSTCluster information to a given
    output stream.  This method provides a convenient mechanism to
    dump the complete MSTCluster information for debugging purposes.
    The clusters in the MSTCluster are displayed in the order in which
    they were created.
*/
extern std::ostream& operator<<(std::ostream&, const MSTCluster&);

#endif
