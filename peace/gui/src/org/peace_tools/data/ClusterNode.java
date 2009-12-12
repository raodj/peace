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

package org.peace_tools.data;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.ProgressMonitor;

/**
 * A class that represents a single cluster in a cluster tree.
 * 
 * This class is a pure data class that is used to encapsulate the information
 * pertaining to a single cluster.  This class is a self-referential structure,
 * in that it contains a list of Cluster objects that represent child clusters of
 * this cluster. This definition permits a Cluster to contain a complete cluster
 * collection a part of it. The Cluster objects are created and used by the 
 * ClusterTree class that represents the top-level class.
 */
public class ClusterNode {
	/** Constructor to create a Cluster.
	 * 
	 * This constructor provides a convenient mechanism to create and 
	 * initialize a Cluster object with necessary information.
	 * 
	 * @param parent The parent cluster for this cluster node (if known).
	 * The parent value is typically set when this node is added as a child
	 * to another cluster node. 
	 * 
	 * @param estNode If this flag is true, that indicates that this
	 * estNode
	 */
	public ClusterNode(ClusterNode parent, boolean estNode, int id) {
		this.parent         = parent;
		this.clusterOrESTID = (estNode ? -(id + 1) : id);
		this.children       = null;
	}

	/**
	 * Add another cluster as the child cluster of this cluster node. This 
	 * method adds the given cluster as a child cluster. In addition it also 
	 * sets up the parent reference in the child to point to this cluster.
	 * 
	 * @param node The cluster node to be added as a direct child of this 
	 * cluster node.
	 */
	public void addChild(ClusterNode node) {
		if (children == null) {
			children = new ArrayList<ClusterNode>();
		}
		node.parent = this;
		children.add(node);
	}
	
	/**
	 * Determine if this cluster node is the root cluster.
	 * 
	 * @note This method is meaningful only after a complete cluster hierarchy
	 * has been built.
	 * 
	 * @return This method returns true if the parent of this cluster is
	 * null, indicating this is a root cluster.
	 */
	public boolean isRoot() { return parent == null; }

	/**
	 * Determine if this cluster is a leaf cluster that has no child clusters
	 * or EST entries.
	 * 
	 * @note This method is meaningful only after a complete cluster hierarchy
	 * has been built.
	 * 
	 * @return This method returns true if this cluster has no children.
	 */
	public boolean isLeaf() { return (children == null); }

	/**
	 * Determine if this node represents a EST entry in the cluster tree.
	 * 
	 * @return This method returns true if this node represents an est
	 * entry in the cluster tree.
	 */
	public boolean isESTNode() { return this.clusterOrESTID < 0; }
	
	/**
	 * Obtain the EST id associated with an EST node.
	 * 
	 * @note The return value from this method is meangiful only if the
	 * isESTNode() method returns true.
	 *  
	 * @return The ID (typically the index of the EST) of the EST 
	 * associated with this cluster node.
	 */
	public int getESTId() { return -(this.clusterOrESTID + 1); }

	/**
	 * Obtain the cluster id associated with an EST node.
	 * 
	 * @note The return value from this method is meangiful only if the
	 * isESTNode() method returns false.
	 *  
	 * @return The ID of this cluster node. 
	 */
	public int getClusterId() { return this.clusterOrESTID; }

	/**
	 * Obtain the child nodes associated with this node.
	 * 
	 * @note This method returns null if the node does not have any
	 * children.
	 * 
	 * @return This method must be used to obtain the child nodes 
	 * associated with a given node. 
	 */
	public ArrayList<ClusterNode> getChildren() { return children; }
	
	/**
	 * Obtain the aggregate classification information for this cluster.
	 * 
	 * This method must be used to obtain the classification information
	 * associated with this cluster. This method returns a HahsMap containing
	 * the classification entries for this cluster. The key values in
	 * the hash map are the index of the DB classifier associated with
	 * a given entry. The value is the number of ESTs in this cluster
	 * that are associated with this EST.
	 * 
	 * @note The return value from this method may be null if a suitable
	 * classifier is not available for this cluster or if a classifier
	 * has not been computed.
	 * 
	 * @return A hash map containing the classification information for
	 * the ESTs. 
	 */
	public HashMap<Integer, Integer> getESTClasses() {
		return estGroups;
	}
	
	/**
	 * This method returns basic information about this cluster node.
	 * 
	 * @return This method returns a simple string representation of
	 * the data stored in this cluster node.
	 */
	@Override
	public String toString() {
		if (!isESTNode()) {
			return "Cluster: " + getClusterId() + " (Nodes: " + 
				((children != null) ? children.size() : 0) + ")";
		}
		return "EST: " + getESTId();
	}

	/**
	 * Determine the largest cluster in this cluster node. 
	 * 
	 * This method recursively searches the cluster hierarchy to 
	 * determine the largest cluster (in terms of ESTs).
	 * 
	 * @return The largest cluster under this cluster node. If the
	 * node is the root, then this return value corresponds to the
	 * largest cluster in the entire cluster file.
	 */
	public int getLargestClusterSize() {
		if (isESTNode()) {
			return 1;
		}
		int estCount      = 0;
		int maxSubCluster = 0;
		for(ClusterNode node: children) {
			if (node.isESTNode()) {
				estCount++;
			} else {
				maxSubCluster = Math.max(maxSubCluster, node.getLargestClusterSize());
			}
		}
		return Math.max(estCount, maxSubCluster);
	}
	
	/**
	 * Method to recursively print the information that is stored in this 
	 * cluster. This method recursively prints the child clusters. This method
	 * also serves as a manual mechanism to validate the data stored in
	 * this cluster. 
	 * 
	 * @param out The output stream to which the data is to be serialized.
	 * @param cluster The cluster to be printed.
	 * @param indent The number of spaces to be used to indent the output.
	 */
	public void print(PrintStream out, ClusterNode cluster, String indent) {
		// First print the EST cluster itself.
		out.println(indent + cluster);
		// Now print the child clusters indented one more step.
		indent = indent + "  ";
		// Print the clusters if we have any.
		if (children != null) {
			for(ClusterNode child: children) {
				child.print(out, child, indent);
			}
		}
	}

	/**
	 * Method to dump the EST data out in FASTA file format. 
	 * 
	 * This method can be used to dump out the ESTs in this cluster to
	 * a file in FASTA compatible format. 
	 * 
	 * @param estList The list of ESTs from which the EST data is to be
	 * obtained for writing.
	 * 
	 * @param os The output stream to which the EST must be written in
	 * a FASTA format.
	 */
	public void write(ESTList estList, PrintStream os) {
		// Print the clusters if we have any.
		if (children != null) {
			for(ClusterNode child: children) {
					child.write(estList, os);
			}
		} else {
			int estIdx = getESTId();
			EST est = estList.getESTs().get(estIdx);
			est.write(os);
		}
	}
	
	/**
	 * Method to compute classification statistics for this cluster node.
	 * 
	 * This method is typically invoked from the ClusterFile.classify()
	 * method to compute classifications for this node. This method 
	 * iterates over all the entries in this cluster and collates 
	 * information about each EST in the cluster. If this cluster has
	 * sub-clusters the the classification is delegated to the sub-clusters
	 * and this node does not really compute any additional classification
	 * information.
	 * 
	 * @param estList The list of ESTs associated with the clusters. This
	 * information is used to compute aggregate classification data for
	 * each cluster.
	 * 
	 * @param pm An optional progress monitor to be used to indicate progress
	 * as the data is computed.
	 */
	protected void classify(ESTList estList, ProgressMonitor pm) {
		// Hash map to track the index of the classifier and the number of entries
		// per classifier in this cluster node.
		estGroups = new HashMap<Integer, Integer>();
		if (children != null) {
			for(int i = 0; (i < children.size()); i++) {
				ClusterNode child = children.get(i);
				if (child.isESTNode() && (estGroups != null)) {
					EST est = estList.getESTs().get(child.getESTId());
					Integer dbClasIdx = new Integer(est.getDBClassifier());
					Integer clasCount = estGroups.get(dbClasIdx);
					if (clasCount == null) {
						// New entry.
						clasCount = new Integer(1);
					} else {
						clasCount = clasCount + 1;
					}
					// Update the group with the latest count
					estGroups.put(dbClasIdx, clasCount);
				} else {
					// Let the child cluster compute its own statistics
					child.classify(estList, null);
					// Stop computing stats as they are not really applicable
					// for this cluster that has sub-clusters in it.
					estGroups = null;
				}
				// Update progress if monitor is not null.
				if (pm != null) {
					pm.setProgress(estList.getESTs().size() + i); 
				}
			}
		}
	}
	
	/** The parent cluster for this MSTcluster. If this MSTcluster is the root 
	 * cluster then it has no parent (that is parent is null).
	 */
	private ClusterNode parent;

	/**
	 * A number that is assigned to this node. This value could either be
	 * the EST index/ID or a cluster ID depending on the sign. Positive
	 * values represent clusterIDs while negative values represent
	 * EST index/ID.
	 */
	private int clusterOrESTID;
	
	/**
	 * This array list contains the list of child clusters and EST nodes
	 * contained by this cluster node. Typically EST nodes will not have
	 * any children. 
	 */
	private ArrayList<ClusterNode> children;
	
	/**
	 * The aggregate classification information associated with this EST.
	 * The information is used to provide additional summary information
	 * about the ESTs in this cluster.
	 */
	HashMap<Integer, Integer> estGroups;
}
