package org.peace_tools.data;

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
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;

import org.netbeans.swing.outline.RowModel;
import org.peace_tools.workspace.FileEntry;

/**
 * The tree table model that provides a facade to display 
 * cluster information in a tree-table.
 * 
 * <p>This class servers as a light weight wrapper to adapt the in-memory 
 * cluster representation to be displayed in a tree-table. The first column
 * in the tree-table displays an hierarchical view of the clusters. The 
 * remaining columns display the EST sequence base pair information. This
 * class serves as the "Model" as in a the Model-View-Controller design
 * pattern.</p> 
 * 
 * <p>The tree-table or Outline model is a generic non-standard GUI component
 * developed by Sun Microsystems as part of the NetBeans IDE. This jar has been
 * obtain from the NetBeans package as a part. The tree-table is a combination 
 * both a JTree and a JTable. The first column in the tree-table is a tree that
 * provides the user with a convenient interface to access and control view of
 * hierarchical information. The remaining columns in the tree-table display 
 * detailed information about each entry in the tree table.</p> 
 */
public class ClusterTreeTableModel implements TreeModel, RowModel {
	/**
	 * The constructor. The constructor merely initializes the necessary
	 * data members in the cluster tree (and its base class).
	 * 
	 * @param clusters The set of clusters to be displayed by this class.
	 * 
	 * @param ests The set of ESTs that contain information about each 
	 * EST in the clusters.
	 * 
	 * @param wsEntry The workspace entry corresponding to the clusters being
	 * adapted by this object.
	 */
	public ClusterTreeTableModel(ClusterFile clusters, ESTList ests, FileEntry wsEntry) {
		// Save references for future use.
		this.clusterFile = clusters;
		this.estList     = ests;
		this.basesPerCol = 100;
		this.wsEntry     = wsEntry;
		// Compute the maximum length of ESTs to determine column count.
		int maxLen = 0;
		for(EST est: estList.getESTs()) {
			maxLen = Math.max(maxLen, est.getSequence().length());
		}
		// Save max len as it it will never change.
		this.maxESTLen = maxLen;
		assert(this.clusterFile != null);
		assert(this.estList != null );
	}

	/**
	 * Obtain the class that describes the data type of a given 
	 * column.
	 * 
	 * This method overrides the API method in RowModel.
	 * 
	 * @param column The zero-based index of the column whose
	 * Class type is to be returned.
	 */
	@Override
	public Class<?> getColumnClass(int column) { 
		return String.class;
	}
	
	/**
	 * This method overrides the interface method in RowModel to 
	 * return the number of columns in the TreeTable. This value is
	 * determined using a combination of maxLen and basesPerCol values.
	 * 
	 *  This method overrides the API method in RowModel.
	 *  
	 * @return The number of columns to be displayed in this table. There
	 * are always a minimum of 1 column.
	 */
	public int getColumnCount() {
		int colCount = (int) Math.ceil(maxESTLen / (double) basesPerCol);
		return colCount;
	}

	/**
	 * This method overrides the interface method in RowModel to 
	 * return the title for the columns. The column titles are dynamically
	 * determined depending on the number of columns displayed in the table.
	 * 
	 * @return The title for the columns. The title for the first column
	 * that contains the tree display is fixed.
	 */	
	public String getColumnName(int column) {
		// Compute title based on the bases per column value.
		return String.format("%d - %d", column * basesPerCol, 
				(column + 1) * basesPerCol);
	}

	/**
	 * This method returns the subset of base pairs to be displayed in a
	 * given column for a given RowModel. This method overrides the 
	 * interface method in TreeTableModel to return the data to be displayed
	 * in a given row. The return value is computed as follows:
	 * 
	 * <ol>
	 * 
	 * <li>If the node is an EST node and column * basesPerCol is less than
	 * EST length, then the bases logically associated with the column are 
	 * returned.</li>
	 * 
	 * <li>Otherwise this method returns an empty string.</li>
	 * 
	 * </ol>
	 * 
	 * @return The data to be displayed for a given node (or row) and a 
	 * given column. 
	 */
	public Object getValueFor(Object node, int column) {
		String      retVal = "";
		ClusterNode entry  = (ClusterNode) node;
		if (entry.isESTNode()) {
			// Get the EST sequence from the ESTList based on EST id
			final EST est    = estList.getESTs().get(entry.getESTId());
			final String seq = est.getSequence();
			final int start  = (column * basesPerCol);
			final int stop   = Math.min(start + basesPerCol, seq.length());
			if (start < seq.length()) {
				retVal = seq.substring(start, stop);
			}
		} else {
			// This is a cluster node. If the column is zero then return the
			// cluster name.
			if (column == 0) {
				retVal = entry.getName(); 
			}
		}
		return retVal;
	}

	/**
	 * Determine if the cell in this column is editable for 
	 * the passed node.
	 * 
	 * This method implements the corresponding method in the
	 * RowTable interface.
	 * 
	 * @param node The node corresponding to the entry in the tree.
	 * @param column The zero-based index of the column to be edited.
	 * @return This method always returns false to indicate that the
	 * data is not editable.
	 */
	@Override
	public boolean isCellEditable(Object node, int column) {
		return false;
	}
	
	/**
	 * Set the value of the object in this column.
	 * 
	 * This method implements the corresponding method in the
	 * RowTable interface. Currently this method does nothing
	 * as the data is not really editable.
	 * 
	 * @param node The node corresponding to the entry in the tree.
	 * @param column The zero-based index of the column to be edited.
	 * @param value The new value to be set for the given column and
	 * node (row).
	 */
	public void setValueFor(Object node, int column, Object value) {
	}
    
	/**
	 * Method to change the order in which top-level clusters are ordered.
	 * 
	 * This method resets the order in which top-level clusters (that is,
	 * clusters that are directly the children of the root) are displayed.
	 * 
	 * @param order The order in which the clusters are to be sorted. The 
	 * following values are valid for this parameter: 0: no sorting,
	 * 1: smaller clusters first, 2: smaller clusters last.
	 * 
	 */
	public void sort(final int order) {
		if (order == 0) {
			// Reset sorting.
			sortedSubClusterIndexs = null;
			// Fire notification
			fireTableStructureChanged();
			return;
		}
		// The comparator object to be used for sorting below.
		Comparator<Integer> comp = new Comparator<Integer>() {
			@Override
			public int compare(Integer param1, Integer param2) {
				ArrayList<ClusterNode> clusters = clusterFile.getRoot().getChildren(); 
				ClusterNode node1   = clusters.get(param1);
				ClusterNode node2   = clusters.get(param2);
				int node1ChildCount = (node1.getChildren() != null) ? node1.getChildren().size() : 1;
				int node2ChildCount = (node2.getChildren() != null) ? node2.getChildren().size() : 1;
				if (order == 1) {
					return node1ChildCount - node2ChildCount;
				}
				return node2ChildCount - node1ChildCount;
			}
		};
		// Create a default list of values.
		ArrayList<ClusterNode> clusters = clusterFile.getRoot().getChildren();
		sortedSubClusterIndexs = new Integer[clusters.size()];
		for(int i = 0; (i < clusters.size()); i++) {
			sortedSubClusterIndexs[i] = new Integer(i);
		}
		// Sort the list based on the sort criteria
		Arrays.<Integer>sort(sortedSubClusterIndexs, comp);
		// Fire notification
		fireTableStructureChanged();
	}
	
	/**
	 * Helper method to broadcast notification to all model listeners.
	 * 
	 * This method is a helper method that is used to broadcast notification
	 * that the tree model has changed to all interested listeners. Typically,
	 * the listeners are views that refresh their display when the data
	 * changes.
	 */
	protected void fireTableStructureChanged() {
		TreeModelEvent tme = new TreeModelEvent(clusterFile.getRoot(),
				new TreePath(clusterFile.getRoot()));
		for(TreeModelListener tml: treeModelListeners) {
			tml.treeStructureChanged(tme);
		}
	}
	
	//------------------------------------------------------------
	
	/**
	 * Adds a listener to be notified when the data associated with
	 * the data set changes.
	 * 
	 * This method implements the API method in TreeModel.
	 * 
	 * @param tml The listener to be added to the list of listeners to
	 * be notified when the data associated with the workspace changes.
	 */
	@Override
	public void addTreeModelListener(TreeModelListener tml) {
		treeModelListeners.add(tml);
	}

	/**
	 * Returns the child cluster or EST node for a given parent node.
	 * 
	 * This method overrides the interface method in the base class to
	 * return the tree node corresponding to the given child. 
	 * This method implements the API method in TreeModel.
	 * 
	 * @return The child cluster or EST node in the given cluster.
	 */
	@Override
	public Object getChild(Object parent, int index) {
		ClusterNode entry  = (ClusterNode) parent;
		if (entry.getChildren() == null) {
			return null;
		}
		if (parent == getRoot() && (this.sortedSubClusterIndexs != null)) {
			// Pay attention to the sorting order in this case.
			index = sortedSubClusterIndexs[index];
		}
		// For nodes whose parent are not root, the order is natural
		// order in which they occur.
		return entry.getChildren().get(index);
	}

	/**
	 * Return the number of child nodes for a given cluster node.
	 * 
	 * This method overrides the interface method in the TreeModel to
	 * return the number of child nodes for a given node.
	 * 
	 * @return The number of child nodes for a given node.
	 */
	@Override
	public int getChildCount(Object parent) {
		ClusterNode entry  = (ClusterNode) parent;
		if (entry.getChildren() != null) {
			return entry.getChildren().size();
		}
		return 0;
	}
	
	/**
	 * Obtain the index of a given child node within its immediate parent.
	 * 
	 * This method overrides the interface method in TreeModel to
	 * return the index of a child node for a given parent node.
	 * 
	 * @return The index of the child node within the parent node.
	 */
	@Override
	public int getIndexOfChild(Object parent, Object child) {
		ClusterNode entry = (ClusterNode) parent;
		if (entry.getChildren() == null) {
			return -1;
		}
		return entry.getChildren().indexOf(child);
	}
	
	/**
	 * Obtain the top-level root node.
	 * 
	 * This method implements the corresponding API method in the
	 * TreeModel. 
	 * 
	 * @return The top-level root cluster node.
	 */
	@Override
	public ClusterNode getRoot() {
		return clusterFile.getRoot();
	}
	
	/**
	 * Determine if a given node in the tree is a leaf node.
	 * 
	 * @param node The node must be a valid node obtained via an
	 * earlier call to the getChild.
	 * 
	 * @return Returns true if the given node is a leaf node in 
	 * the tree.
	 */
	@Override
	public boolean isLeaf(Object node) {
		ClusterNode clsNode = (ClusterNode) node;
		return clsNode.isLeaf();
	}
	
	/**
	 * Removes a listener from the list of listeners maintained by
	 * this class.
	 * 
	 * @param tml The listener to be removed from the list of listeners
	 * to receive updates.
	 */
	@Override
	public void removeTreeModelListener(TreeModelListener tml) {
		treeModelListeners.remove(tml);
	}

	/**
	 * Messaged when the user has altered the value for the 
	 * item identified by path to newValue.
	 *
	 * This method is currently not implemented because the cluster
	 * data cannot be edited.
	 */
	public void valueForPathChanged(TreePath path, Object newValue) {
		// Nothing to be done for now.
	}
	
	/**
	 * Obtain the list of ESTs associated with this model.
	 * 
	 * @return The EST list set for use by this model.
	 */
	public ESTList getESTList() { return estList; }

	/**
	 * Obtain the cluster file from where cluster data was obtained.
	 * 
	 * @return The cluster file from where the clustering information
	 * was obtained.
	 */
	public ClusterFile getClusterFile() { return clusterFile; }
	
	/**
	 * Obtain the actual workspace entry whose data is contained in this
	 * model.
	 * 
	 * A handy reference to the workspace entry from which the data for this
	 * cluster table model was actually obtained. This information can be
	 * used by "view" classes to create additional views as needed.
	 * 
	 * @return The reference to the MSTClusterData workspace entry whose
	 * data is "modeled" by this class.
	 */
	public FileEntry getWsEntry() { return wsEntry; }
	
	//-------------------------------------------------------
	

	/**
	 * Reference to the cluster file that contains the cluster data to
	 * be adapted and interfaced by this model. This value is set in the
	 * constructor and is never changed during the life time of this class.
	 */
	private final ClusterFile clusterFile;
	
	/**
	 * A handy reference to the workspace entry from which the data for this
	 * cluster table model was actually obtained. This information can be
	 * used by "view" classes to create additional views as needed.
	 */
	private final FileEntry wsEntry;
	
	/**
	 * Reference to the EST file that contains the cluster data to
	 * be adapted and interfaced by this model. This value is set in the
	 * constructor and is never changed during the life time of this class.
	 */
	private final ESTList estList;
	
	/**
	 * The maximum length of an EST sequence to be adapted by this model.
	 * This value along with the basesPerCol determines the total number of 
	 * columns that are logically represented by this model. 
	 */
	private final int maxESTLen;
	
	/**
	 * The number of base pairs per column. This value along with the
	 * maxESTLen determines the total number of columns that are logically
	 * represented by this model. 
	 */
	private int basesPerCol;
	
	/**
	 * This array contains a sorted list of indexes of fragments if a sorting
	 * scheme has been applied to the data set. If this array is not null then
	 * this table model returns entries using the order specified in this
	 * array.  The entries in this array are indexes into the original set 
	 *  of clusters in the root. This list changes depending on the order 
	 *  of sorting requested by the user.
	 */
	private Integer[] sortedSubClusterIndexs = null;
	
    /**
     * The list of tree mode listeners that were added to this tree
     * model. This list is used to notify listeners that the tree
     * has changed.
     */
    private ArrayList<TreeModelListener> treeModelListeners =
	        new ArrayList<TreeModelListener>();
}
