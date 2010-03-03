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
import java.util.Collections;
import java.util.Comparator;

import javax.swing.ComboBoxModel;
import javax.swing.ListModel;
import javax.swing.event.ListDataEvent;
import javax.swing.event.ListDataListener;

import org.peace_tools.workspace.MSTClusterData;

/**
 * A custom data model to provide pre-alignment information using MST and
 * clustering data.
 * 
 * <p>This class servers as a light weight wrapper to adapt the in-memory 
 * MST and clustering representation to be displayed in a custom view. The 
 * term "model" is being used in the context of the Model-View-Controller 
 * (MVC) pattern. The custom view provides a high-level summary information
 * about the type of alignment that is to be expected from the clustering 
 * data. The MST and clustering data are used in the following manner:
 * 
 * <ul>
 * 
 * <li><p>The MST data contains information as to where the maximum overlap
 * between a parent-and-child was detected. This information is used to
 * layout the fragments in a FASTA file in a tabular form with related
 * fragments organized next to each other illustrating the overlap.</p> 
 * 
 * <p>Note that not all MST files contain the relative location information.
 * This information is vital for this data model to work. Consequently, when
 * MST data does not have the positional information, this model generates 
 * an error.</p></li> 
 * 
 * <li>The clustering information is essentially used to color code the 
 * fragments based on the clusters they belong to. In addition, the clustering
 * information provides the user with a convenient mechanism to selectively
 * highlight clusters and study their organization.</li>
 * 
 * </ul>
 * 
 * This data model combines and provides access to the following two distinct 
 * sets of information:
 * 
 * <ol>
 * 
 * <li>First it provides a list of clusters and associated coloring information
 * for each cluster. This information is exposed via the ComboBoxModel interface.</li>
 * 
 * <li>It exposes the data required to draw various EST entries in a 2-dimensional
 * matrix kind-of a layout. </li>
 *  
 * </ol>
 */
public class OverlapModel implements ComboBoxModel {
	/**
	 * The constructor. The constructor merely initializes the necessary
	 * data members in this data model.
	 * 
	 * @param clusters The set of clusters to be displayed by this class.
	 * 
	 * @param ests The set of ESTs that contain information about each 
	 * EST in the clusters.
	 */
	public OverlapModel(ClusterFile clusters, ESTList ests, MSTClusterData wsEntry) {
		this.estList         = ests;
		this.clusterFile     = clusters;
		this.selectedCluster = -1;
	}

	/**
	 * The constructor. The constructor merely initializes the necessary
	 * data members in this data model.
	 * 
	 * @param clusters The set of clusters to be displayed by this class.
	 * 
	 * @param ests The set of ESTs that contain information about each 
	 * EST in the clusters.
	 */
	public static OverlapModel create(ClusterFile clusters, ESTList ests, MST mst) {
		if (!mst.hasAlignmentInfo()) {
			// Can't use this MST as there is no alignment 
			return null;
		}
		// Create the classs with necessary information.
		OverlapModel pam    = new OverlapModel(clusters, ests, null);
		// Create the list of clusters to operate on.
		pam.buildClusterList(clusters);
		int[] clusterIDList = new int[ests.getESTs().size()];
		OverlapModel.generateClusterLookupTable(-1, clusters.getRoot(), clusterIDList);
		int leftMostPos          = pam.getFarthestPos(mst.getRoot());
		ArrayList<Long> rowUsage = new ArrayList<Long>();
		pam.generateESTEntries(mst.getRoot(), -leftMostPos, clusterIDList, -1, -1, rowUsage);
		// Everything went well. Return the newly created model for further use
		return pam;
	}

	
	/**
	 * Obtain the currently displayed cluster in the combo box.
	 * 
	 * This method can be used to obtain the current cluster being displayed
	 * (or to be displayed) in a combo box. This method implements the
	 * corresponding interface method defined in the ComboBoxModel.
	 * 
	 * @return This method returns the currently selected item in the 
	 * list.
	 */
	@Override
	public ClusterNode getSelectedItem() {
		if ((selectedCluster >= 0) && (selectedCluster < clusterList.size())) {
			return clusterList.get(selectedCluster);
		}
		// Invalid index.
		return null;
	}


	/**
	 * Set the selected (item being displayed) item in a combo box. 
	 * 
	 * As per the API requirement this method should notifies all registered
	 * ListDataListeners that the contents have changed.
	 * 
	 *  @param anItem The item to be displayed in the combo box.
	 * 
	 */
	@Override
	public void setSelectedItem(Object anItem) {
		selectedCluster = clusterList.indexOf(anItem);
		// Report the fact that things have changed.
		fireDataChanged();
	}
	
	/**
	 * Adds a listener to be notified when the data associated with
	 * the model changes. This method is part of the ListModel
	 * interface class.
	 * 
	 * @param ldl The listener to be added to the list of listeners to
	 * be notified when the data associated with the workspace changes.
	 * 
	 * @see ListModel
	 */
	@Override
	public void addListDataListener(ListDataListener ldl) {
		listListeners.add(ldl);
	}

	/**
	 * Removes a listener from the list of listeners maintained by
	 * this class. This method is part of the ListModel
	 * interface class.
	 * 
	 * @param ldl The listener to be removed from the list of listeners
	 * to receive updates.
	 * 
	 * @see ListModel
	 */
	@Override
	public void removeListDataListener(ListDataListener ldl) {
		 listListeners.remove(ldl);
	}

	/**
	 * Obtain the cluster node in the combo box list at a given index.
	 * 
	 * This method implements the interface method in ListModel to
	 * return the appropriate cluster entry at a given index. Note that
	 * the actual cluster node returned may vary depending on how the
	 * list is sorted.
	 * 
	 * @param index The index of the entry whose element is to be returned
	 * by this method. This value must be in the range 0 &le; index &lt; getSize()
	 * 
	 * @return The element at the requested index. If the index is not valid then
	 * this method returns null.
	 */
	@Override
	public ClusterNode getElementAt(int index) {
		if ((index >= 0) && (index < clusterList.size())) {
			return clusterList.get(index);
		}
		// Invalid index.
		return null;
	}

	/**
	 * Obtain the number of entries to be displayed in the combo box.
	 * 
	 * This method implements the corresponding interface method in the ListModel.
	 * It returns the number of clusters involved in this data model. Note that
	 * clusters with names (that currently represent dummy clusters) are not
	 * included in this list.
	 * 
	 * @return The number of actual clusters to be displayed in this model.
	 */
	@Override
	public int getSize() {
		return clusterList.size();
	}
	
	/**
	 * Obtain the list of ESTs associated with this model.
	 * 
	 * @return The EST list set for use by this model.
	 */
	public ESTList getESTList() { return estList; }
	
	/**
	 * Obtain the cluster file information from where this 
	 * pre-alignment model was built.
	 * 
	 * @return The cluster file (containing cluster information and
	 * other meta data) set for use by this model.
	 */
	public ClusterFile getClusterFile() { return this.clusterFile; }

	/**
	 * Obtain the last column where a nucleotide is present.
	 * 
	 * This method can be used to determine the logical width of this model.
	 * This value is typically used to determine the overall viewing area
	 * required to display the model.
	 * 
	 * @return The last column where a nucleotide is present. This value is
	 * at least 1.
	 */
	public int getMaxCol() { return maxCol; }
	
	/**
	 * Determine the maximum number of rows in this overlap model.
	 * 
	 * This method must be used to determine the maximum number of rows in
	 * this overlap model.
	 * 
	 * @return The maximum number of rows in this overlap model.
	 */
	public int getMaxRow() {
		return this.estEntryTable.size();
	}
	
	/**
	 * Obtain the list of EST entries on a given row in the model.
	 * 
	 * This method may be used to obtain the list of fragment entries on
	 * a specific row of the model.
	 * 
	 * @param row The row for which the list of fragments are to be returned.
	 * 
	 * @return An array list containing the list of fragments on a given row. If
	 * the row was invalid, then this method return null.
	 */
	public ArrayList<ESTEntry> getRow(final int row) {
		if ((row >= 0) && (row < estEntryTable.size())) {
			return estEntryTable.get(row);
		}
		return null;
	}
	
	//-------------------------------------------------------
	
	/**
	 * Helper method to notify listeners whenever data changes.
	 * 
	 * This method is used by other methods in this class to notify all
	 * list listeners whenever the data associated with the combo model
	 * changes. This method currently reports that the contents has
	 * changed to all registered listeners.
	 */
	protected void fireDataChanged() {
		ListDataEvent event = new ListDataEvent(this, ListDataEvent.CONTENTS_CHANGED, 
				0, clusterList.size()); 
		for(ListDataListener listener: listListeners) {
			listener.contentsChanged(event);
		}
	}
	
	/**
	 * Method to change the order in which clusters are ordered.
	 * 
	 * This method resets the order in which clusters are delivered to
	 * the view using this model.
	 * 
	 * @param order The order in which the clusters are to be sorted. The 
	 * following values are valid for this parameter: 0: no sorting,
	 * 1: smaller clusters first, 2: smaller clusters last.
	 * 
	 */
	public void sort(final int order) {
		// The comparator object to be used for sorting below.
		Comparator<ClusterNode> comp = new Comparator<ClusterNode>() {
			@Override
			public int compare(ClusterNode node1, ClusterNode node2) {
				int node1ChildCount = (node1.getChildren() != null) ? node1.getChildren().size() : 1;
				int node2ChildCount = (node2.getChildren() != null) ? node2.getChildren().size() : 1;
				int retVal          = (node1.getClusterId() - node2.getClusterId());
				if (order != 0) {
					retVal = (order == 1) ? (node1ChildCount - node2ChildCount) :
						(node2ChildCount - node1ChildCount);
				} 
				return retVal;
			}
		};
		// Sort the list of clusters based on the sort criteria
		Collections.sort(clusterList, comp);
		// Fire notification
		fireDataChanged();
	}
	
	/**
	 * Helper method to determine how far-left the pre-alignment extends given
	 * the root is at logical position zero.
	 * 
	 * This method recursively calls itself to determine the left-most 
	 * fragment position.  Essentially, this method first determines the
	 * farthest positions of each child node (given the parent node) and 
	 * returns the minimum of the position values.  The recursion terminates
	 * each time a leaf node in the MST is encountered.
	 * 
	 * @param parent The parent node whose farthest child position is to be
	 * determined. This method is typically invoked with the root of the MST 
	 * tree.
	 * 
	 * @return This method returns the base-pair position of the lest-most
	 * fragment in the MST
	 */
	private int getFarthestPos(final MSTNode parent) {
		int farthestPos = 0;
		if ((parent == null) || (parent.getChildren() == null)) {
			// No further checks needed.
			return 0;
		}
	    // Find the minimum of the farthestPos amongst all the child 
		// nodes in a recursive manner.
		for(MSTNode childNode: parent.getChildren()) {
			// Obtain child's (sub-tree's) estimate of farthest position.
			final int childFarthestPos = getFarthestPos(childNode);
			// Compute the column where the child is going to be placed.
			final int sign = (Math.signum(childNode.getAlignmentMetric()) >= 0 ? 1 : -1);
			final int childPos = (int) (childNode.getAlignmentMetric() + 
				(childNode.getMetric() * sign)) + childFarthestPos; 
			// Track the minimum we have computed thus far
			farthestPos = Math.min(farthestPos, childPos);
		}
		// Return the minimum estimate.
		return farthestPos;
	}
	
	/**
	 * Helper method to populate local list of clusters from a cluster file.
	 * 
	 * This method is a helper method that is used to process a given 
	 * ClusterFile and add non-dummy cluster nodes to the local list of
	 * clusters maintained by this class. This method is used only once
	 * while the PreAlignmentModel is being populated.
	 * 
	 * @param clusterFile The cluster file to be processed by this method to
	 * extract the list of clusters to be managed by this class.
	 */
	private void buildClusterList(ClusterFile clusterFile) {
		// Obtain the root node from the clusterFile.
		final ClusterNode root = clusterFile.getRoot();
		// Clear and reserve space in our local cluster list.
		clusterList = new ArrayList<ClusterNode>();
		clusterList.ensureCapacity(root.getChildren().size());
		// Add non-dummy clusters from root node to our local list
		for(ClusterNode child: root.getChildren()) {
			if (child.isESTNode() || !child.getName().isEmpty()) {
				// This entry is either a EST node or a dummy cluster.
				// So ignore this entry.
				continue;
			}
			// Found a valid cluster. Add it to the local list
			clusterList.add(child);
		}
	}
	
	/**
	 * Add an entry to the specified row in the estEntryTable.
	 * 
	 * This method is a helper method that is invoked from the generateESTEntries
	 * method to add a new entry to the estEntryTable. This method also
	 * tracks and updates {@link OverlapModel#maxCol} instance variable
	 * to reflect the farthest nucleotide in the table (across all rows).
	 * 
	 * @param row The row in the estEntryTable to which the given entry is
	 * to be added.
	 * 
	 * @param entry The newly created entry to be added to the estEntryTable
	 * at the specified row.
	 */
	private void addToESTEntryTable(final int row, ESTEntry entry) {
		// First ensure we have suffient number of (emtpy) rows in the table
		while (estEntryTable.size() <= row) {
			estEntryTable.add(new ArrayList<ESTEntry>());
		}
		// Now add the entry at the appropriate row.
		estEntryTable.get(row).add(entry);
		// Update maxCol as agreed by the API.
		final int estLen = entry.getEST().getSequence().length();
		maxCol = Math.max(maxCol, entry.getColumn() + estLen);
	}
	
    /** Recursive method to generate ESTEntry objects with alignment information.
     * 
     * This method is invoked to perform the actual core alignment task using the
     * input fragments from a FASTA file along with the MST (with alignment 
     * information).  This method uses the startPos value provided by the caller
     * to recursively position fragments and generate ESTEntry objects with 
     * alignment information for each fragment.
     * 
     * @param parent The MST node (and its children) that must be aligned
     * by this method. This method is initially invoked on the root of the MST.
     * 
     * @param startPos The starting base pair position relative to which the
     * parent and its child nodes are to be positioned.
     * 
     * @param clusterIDList The lookup table that has been populated by the
     * {@link OverlapModel#generateClusterLookupTable(int, ClusterNode, int[])}
     * method. This method uses the look up table to determine the cluster ID when
     * generating ESTEntry objects for each fragment.
     * 
     * @param prevClusterID This parameter is used to track the ID of the 
     * previous cluster to which an entry was added. When moving from one cluster
     * to another, the lastRow value is reset to 0.
     * 
     * @param lastRow The last row in the rowUsage list where the previous entry was
     * placed.
     * 
     * @param rowUsage The array of row usage information that is used to determine
     * when a row is available to accommodate a given entry.
     */    
	private void generateESTEntries(final MSTNode parent, final int startPos, 
			final int clusterIDList[], int prevClusterID, 
			int lastRow, ArrayList<Long> rowUsage) {
		// First generate a suitable entry for the parent node.
		// Look up the clusterID for the parent 
		final int estIdx    = parent.getESTIndex();
		final int clusterID = clusterIDList[estIdx];
		if (clusterID != -1) {
			final EST est       = estList.getESTs().get(estIdx);
			// Create the entry for the parent node.
			ESTEntry parEntry   = new ESTEntry(est, clusterID, startPos);
			// Add it to our list of EST entries
			estEntryList.add(parEntry);
			// Next place it in an appropriate row that is adjacent to the previous one.
			if (prevClusterID != clusterID) {
				prevClusterID = clusterID;
				lastRow       = 0;
			}
			// Get the row were the entry should be placed.
			lastRow = getRow(startPos, startPos + est.getSequence().length(), 
							 lastRow - 1, rowUsage);
			// Add entry to the specified row.
			addToESTEntryTable(lastRow, parEntry);
		}
		
		if (parent.getChildren() == null) {
			// This node does not have any child nodes.
			// Nothing further to be done. 
			return;
		}
		// Now let child nodes generate entries for themselves.
		for(MSTNode child: parent.getChildren()) {
			// First compute the logical column position for the child node
			// with respect to the parent position 
			final int sign = (Math.signum(child.getAlignmentMetric()) >= 0 ? 1 : -1);
			final int childPos = startPos + child.getAlignmentMetric() +
				(int) (child.getMetric() * sign);
			// Let the child create its own entry.
			generateESTEntries(child, childPos, clusterIDList, prevClusterID, lastRow, rowUsage);
		}
	}

	/**
	 * Helper method to build a table (simple array) that provides clusterID for each
	 * fragment.
	 * 
	 * This method is a recursive helper method that is used to build a look-up table
	 * to map a fragment to its cluster ID. The look-up table is simply an integer array.
	 * The index into this array is the index of the fragment whose cluster ID is to
	 * be determined. The array contains the cluster ID for each fragment read from a
	 * given FASTA file.
	 *  
	 * @param parentClusterID The logical cluster ID of the parent cluster to which the
	 * given child node logically belongs.
	 * 
	 * @param node The child node in the cluster hierarchy whose leaf fragment entries
	 * are to be cataloged in the lookup table.
	 * 
	 * @param clusterIDList The lookup table being populated by this method. The array must
	 * be exactly the same size as the number of fragments in the FASTA file.
	 * 
	 * @see OverlapModel#generateESTEntries(MSTNode, int)
	 */
	private static void generateClusterLookupTable(int parentClusterID, ClusterNode node, 
			int clusterIDList[]) {
		if (node.isESTNode()) {
			// This is a leaf node. Update clusterIDList look up table.
			int estIdx = node.getESTId();
			assert(estIdx >= 0);
			assert(estIdx < clusterIDList.length);
			clusterIDList[estIdx] = parentClusterID; 
		} else {
			// This is a cluster entry that has leaf nodes in it. So
			// let them recursively process their child nodes. For dummy
			// clusters flag their clusterID as -1.
			final int clusterID = (node.getName().isEmpty() ? node.getClusterId() : -1);
			if ((node.getChildren() != null) && (!node.getChildren().isEmpty())) {
				for(ClusterNode child: node.getChildren()) {
					generateClusterLookupTable(clusterID, child, clusterIDList);
				}
			}
		}
	}

	/**
	 * A simple helper method to combine two integers into a long.
	 * 
	 * This is a helper method that is used in the getRow method to 
	 * combine two integers into a long value.
	 * 
	 * @param hiInt The high 32-bit value to be stored in the long.
	 * @param lowInt The low 32-bit value to be stored in the long.
	 * @return The long that contains the hiInt and lowInt values combined
	 * into one.
	 */
	private static long toLong(int hiInt, int lowInt) {
		long value = hiInt;
		value    <<= 32;
		value     |= lowInt;
		return value;
	}
	
	/**
	 * This is a helper method that is used to determine the row on which a 
	 * fragment can be accommodated.
	 * 
	 * @param start The starting column where the fragment is to be placed.
	 * @param end The ending column where the fragment will end.
	 * @param startRow The starting row from where the search is to be done.
	 * This value must be 0.
	 * @param rowUsage The row 
	 * @return The row where the fragment can be accommodated.
	 */
	private static int getRow(final int start, final int end, final int startRow, 
			ArrayList<Long> rowUsage) {
		assert ( start < end );
		// These variables are used below the for-loop
		int row = Math.max(0, startRow); // Ensure row is no smaller than 0
		int rowStart = 0, rowEnd = 0;
		// Search existing rowUsage to see if a free row is available
		for(; (row < rowUsage.size()); row++) {
			// Get the starting and ending column of the row usage
			// that has been packed into a Long.
			rowStart = (int) rowUsage.get(row).intValue();
			rowEnd   = (int) (rowUsage.get(row).longValue() >> 32);
			// Check if the existing row can be reused
			if ((rowEnd < start) || (rowStart > end)) {
				// This row can be reused
				break;
			}
		}
		// Check and reuse existing row or add new entry
		if (row == rowUsage.size()) {
			// An existing row cannot be reused. Add a new row.
			rowUsage.add(new Long(toLong(end, start)));
		} else {
			// Update the usage column for the row
			rowStart = Math.min(start, rowStart);
			rowEnd   = Math.max(rowEnd, end);
			// Save information for later use
			rowUsage.set(row, new Long(toLong(rowEnd, rowStart)));
		}
		// Return the row to use 
		return row;
	}
	
	//-------------------------------------------------------

	/**
	 * Reference to the EST file that contains the fragment data that is to
	 * be adapted and presented by this model. This value is set in the
	 * constructor and is never changed during the life time of this class.
	 */
	private final ESTList estList;
	
	/**
	 * Reference to the cluster file that contains the cluster information
	 * that was used to generate this pre-alignment model. This value is set
	 * in the constructor and is never changed during the life time of this 
	 * object.
	 */
	private final ClusterFile clusterFile;
	
	/**
	 * A table (A sparse 2-dimensional array) of ESTEntries.
	 * 
	 * This array list is organized as a 2-dimensional array of ESTEntries.
	 * Each entry in the estEntryTable corresponds to a set of entries to be
	 * displayed as a single row. Accordingly, each row consists of a list 
	 * of ESTEntries. The array is populated by the 
	 * {@link OverlapModel#generateESTEntries(MSTNode, int, int[], int, int, ArrayList)}
	 * method.
	 * 
	 * <b>Note</b> The ESTEntry objects in each row of this table are sorted
	 * based on the starting column of each entry. This information is used by
	 * the view to optimize rendering.
	 */
	private ArrayList<ArrayList<ESTEntry>> estEntryTable = 
		new ArrayList<ArrayList<ESTEntry>>();
	
	/**
	 * The widest row in the estEntryTable. This instance variable is
	 * used to track the widest row in the estEntryTable. This value indicates
	 * the highest column where a nucleotide from some fragment needs to be
	 * rendered. This value is computed and updated each time a new
	 * entry is added to the estEntryTable.
	 */
	private int maxCol = 0;
	
	/**
	 * The list of EST entry objects encapsulated by this data model.
	 * This list of objects is generated by the generateESTEntries method
	 * when an instance of this data model is created.
	 * 
	 * @see OverlapModel#generateESTEntries(MSTNode, int)
	 */
	private ArrayList<ESTEntry> estEntryList = new ArrayList<ESTEntry>();
	
	/**
	 * This list contains the list of clusters (that may be sorted) being
	 * displayed to the user. This list is populated initially when this
	 * data model is created. After that the order of entries in this list
	 * may change whenever the list is (re)sorted.
	 */
	private ArrayList<ClusterNode> clusterList;
	
    /**
     * The set of list model listeners that were added to this model. 
     * This list is used to notify listeners whenever the data in the
     * list model has changed. The data changes whenever the user chooses
     * to resort the list of clusters based on their cluster size.
     */
    private ArrayList<ListDataListener> listListeners =
	        new ArrayList<ListDataListener>();
    
    /**
     * Instance variable to track the currently selected/visible 
     * cluster entry in a combo box. This instance variable is used to
     * implement the functionality associated with the two methods
     * defined in the ComboBoxModel interface class.
     */
    private int selectedCluster;
}
