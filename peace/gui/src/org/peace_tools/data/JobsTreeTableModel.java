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

import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;

import org.netbeans.swing.outline.RowModel;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobList;
import org.peace_tools.workspace.JobList.JobOrListWrapper;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;
import org.peace_tools.workspace.WorkspaceEvent;
import org.peace_tools.workspace.WorkspaceEvent.EntryType;

/**
 * The Job tree table model that provides a facade to display 
 * job hierarchy information in a tree-table.
 * 
 * <p>This class servers as a light weight wrapper to adapt the in-memory 
 * Job hierarchy representation to be displayed in a tree-table. The 
 * first column in the tree-table displays an hierarchical view of 
 * the jobs. The remaining columns display information about the job.</p> 
 * 
 * <p>The tree-table or Outline model is a generic non-standard GUI component
 * developed by Sun Microsystems as part of the NetBeans IDE. This jar has been
 * obtain from the NetBeans package as a part. The tree-table is a combination 
 * both a JTree and a JTable. The first column in the tree-table is a tree that
 * provides the user with a convenient interface to access and control view of
 * hierarchical information. The remaining columns in the tree-table display 
 * detailed information about each entry in the tree table.</p> 
 */
public class JobsTreeTableModel extends JobListTableModel  
implements TreeModel, RowModel {
	/**
	 * The root of the job hierarchy that contains the jobs and job
	 * sublists in it. This value is set when the view associated
	 * with this model is created.
	 */
	private JobOrListWrapper root;

	/**
	 * The list of tree model listeners that were added to this 
	 * model. This list is used to notify listeners that the tree
	 * has changed whenever entries are added, removed, or
	 * reorganized.
	 */
	private ArrayList<TreeModelListener> treeModelListeners =
		new ArrayList<TreeModelListener>();

	/**
	 * The only constructor for this class.
	 * 
	 * The constructor merely initializes the instance variables
	 * to their default initial values.
	 * 
	 * @param root The top-level job list from where this model must
	 * obtain the necessary information.
	 */
	public JobsTreeTableModel(final JobList root) {
		super(root);
		this.root = root.new JobOrListWrapper(root);
	}

	/**
	 * Obtain the class that describes the data type of a given 
	 * column.
	 * 
	 * This method overrides the corresponding method in the base
	 * class to account for the fact that the zero-index column in 
	 * a tree-table view is actually the first column in the table
	 * view as the job ID is displayed in the tree (which is
	 * not considered as the first column)
	 * 
	 * @param column The zero-based index of the column whose
	 * Class type is to be returned.
	 */
	@Override
	public Class<?> getColumnClass(int column) {
		return super.getColumnClass(column + 1);
	}

	/**
	 * This method overrides the interface method in RowModel to 
	 * return the number of columns in the TreeTable.
	 * 
	 * This method implements the corresponding API method in RowModel.
	 * It also overrides the corresponding method in the parent
	 * class to account for the minor difference between the parent's
	 * table model and this tree-table model.
	 *  
	 * @return The number of columns to be displayed in a tree-table. There
	 * are always 4 columns to be displayed.
	 */
	@Override
	public int getColumnCount() {
		return 4;
	}

	/**
	 * Override the default names that are set for the columns 
	 * displayed in the job table.
	 * 
	 * This method uses the base class method with column index
	 * suitably adjusted to account for the differences between
	 * a traditional table-model and this tree-table model.
	 * 
	 * @param col The zero-based column index whose title is to 
	 * be returned.
	 * 
	 * @return The title associated with the column.
	 */
	@Override
	public String getColumnName(int col) {
		return super.getColumnName(col + 1);
	}

	/**
	 * This method returns different information about a Job or JobList
	 * thereby providing a convenient interface for views.
	 * 
	 * <p><b>NOTE</b>: There is a difference in the values returned
	 * by this method and the {@link #getValueAt(int, int)} method
	 * for different columns.</p> 
	 * 
	 * This method returns one of the following types of values depending
	 * on the type of object and the column.
	 * 
	 * <ul>
	 * 
	 * <li>If the node corresponds to a Job node then this method 
	 * returns the following information (consistent with the column
	 * titles returned by {@link #getColumnName(int)}):
	 * 
	 * <ol>
	 * 	<li>If the column number is negative (a special case that normally does
	 * not occur in a GUI but a convenience for views to use) then this
	 * method returns the wrapper entry itself. This value is currently
	 * used by {@link org.peace_tools.views.JobsTreeTableView} to display
	 * a progress bar for running jobs</li>
	 * 
	 * <li>For other column numbers, this method returns the information
	 * consistent with the column titles returned by {@link #getColumnName(int)}.
	 * </ol>
	 * 
	 * </li>
	 * 
	 * <li>For nodes that correspond to a JobList node and for column values
	 * less than zero this method returns the node itself as convenience.
	 * For all other column numbers this method returns an empty string.</li>
	 * 
	 * <li>Otherwise this method returns an empty string.</li>
	 * 
	 * </ul>
	 * 
	 * @return The data to be displayed for a given node (or row) and a 
	 * given column. 
	 */
	@Override
	public Object getValueFor(Object node, int column) {
		String      retVal = "";
		JobOrListWrapper wrapper  = (JobOrListWrapper) node;
		if (column < 0) {
			return wrapper;
		}
		if (!wrapper.isSubList()) {
			Workspace ws = Workspace.get();
			Job job      = wrapper.getJob();
			Server srvr  = ws.getServerList().getServer(job.getServerID());
			switch (column) {
			case 0: // return the current status.
				return job.getStatus();
			case 1: // return monitor status.
				if (job.getMonitor() != null) {
					return "Running";
				} else if (!job.isDone() && !job.isWaiting()) {
					return "Needed";
				}
				return "Not Needed";
			case 2: // return the server on which job is running.
				return (srvr != null) ? srvr.getName() : "<n/a>";
			case 3: // return number of CPUs 
				return "" + (job.getCPUsPerNode() * job.getNodes());
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
		// No, the cell is not editable
		return false;
	}

	/**
	 * Set the value of the object in this column.
	 * 
	 * This method implements the corresponding method in the
	 * RowTable interface. Currently this method does nothing
	 * as the data is read-only.
	 * 
	 * @param node The node corresponding to the entry in the tree.
	 * @param column The zero-based index of the column to be edited.
	 * @param value The new value to be set for the given column and
	 * node (row).
	 */
	@Override
	public void setValueFor(Object node, int column, Object value) {
		// Nothing to be done here.
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
	public JobOrListWrapper getRoot() {
		return root;
	}

	/**
	 * Returns a JobOrListWrapper child node for the given parent
	 * node.
	 * 
	 * This method overrides the interface method in the base class to
	 * return the tree node corresponding to the given child. 
	 * This method implements the API method in TreeModel.
	 * 
	 * @return The child node node in the given job list.
	 */
	@Override
	public JobOrListWrapper getChild(Object parent, int index) {
		JobOrListWrapper entry  = (JobOrListWrapper) parent;
		if (!entry.isSubList()) {
			// Job entries do not have any children
			return null;
		}
		return entry.getSubList().getJobs().get(index);
	}

	/**
	 * Return the number of child nodes for a given node.
	 * 
	 * This method overrides the interface method in the TreeModel to
	 * return the number of child nodes for a given node.
	 * 
	 * @return The number of child nodes for a given node.
	 */
	@Override
	public int getChildCount(Object parent) {
		JobOrListWrapper entry  = (JobOrListWrapper) parent;
		if (entry.isSubList()) {
			return entry.getSubList().getJobs().size();
		}
		// Job nodes do not have any children and its child count is zero.
		return 0;
	}

	/**
	 * Determine if a given node in the tree is a leaf node.
	 * 
	 * All job nodes are leaf nodes while all JobList are 
	 * non-leaf nodes.
	 * 
	 * @param node The node must be a valid node obtained via an
	 * earlier call to the getChild.
	 * 
	 * @return Returns true if the given node is a leaf node in 
	 * the tree.
	 */
	@Override
	public boolean isLeaf(Object node) {
		JobOrListWrapper entry  = (JobOrListWrapper) node;
		return !entry.isSubList();
	}

	/**
	 * Messaged when the user has altered the value for the 
	 * item identified by path to newValue.
	 *
	 * This method is currently not implemented because the job
	 * data cannot be edited.
	 */
	@Override
	public void valueForPathChanged(TreePath path, Object newValue) {
		// Nothing to be done for now
	}

	/**
	 * Obtain the index of a given child node within its 
	 * immediate parent.
	 * 
	 * This method overrides the interface method in TreeModel to
	 * return the index of a child node for a given parent node.
	 * 
	 * @return The index of the child node within the parent node.
	 */
	@Override
	public int getIndexOfChild(Object parent, Object child) {
		JobOrListWrapper entry = (JobOrListWrapper) parent;
		if (!entry.isSubList()) {
			return -1;
		}
		return entry.getSubList().getJobs().indexOf(child);
	}

	/**
	 * Helper method recursively search for a given object starting 
	 * with a given parent node in the model tree.
	 * 
	 * This is a helper method that is invoked from the 
	 * {@link #workspaceChanged(WorkspaceEvent)} method to determine
	 * the full path to the tree node that has changed. This method
	 * recursively calls itself to locate the path for the given
	 * entry. Once the entry is found, the recursion unwinds. As
	 * the methods unwind entries are added to the outgoing
	 * path array list. 
	 * 
	 * @param parent The node from where this class should recursively
	 * search for the given entry.
	 * 
	 * @param entry The entry to search for. This entry can either be a
	 * JobOrListWrapper object or a Job or JobList entry encapsulated
	 * within a given wrapper.
	 * 
	 * @param path The outgoing array list to which the path entries are
	 * to be added. This parameter cannot be null. Entries are added
	 * to this array list from leaf to root.
	 * 
	 * @return This method returns true if the entry was found. Otherwise
	 * it returns false.
	 */
	private boolean addToPath(JobOrListWrapper parent, Object entry,
			ArrayList<JobOrListWrapper> path) {
		if ((parent == entry) || 
			(parent.isSubList()  && (parent.getSubList() == entry)) || 
			(!parent.isSubList() && (parent.getJob()     == entry))) {
			path.add(parent);
			return true;
		}
		// If the parent is a sublist then recursively search within sulist.
		if (parent.isSubList()) {
			for(JobOrListWrapper jlw: parent.getSubList().getJobs()) {
				if (addToPath(jlw, entry, path)) {
					path.add(parent);
					return true;
				}
			}
		} 
		return false;
	}

	/**
	 * Broadcast changes to the workspace-wide unique job list to 
	 * all registered listeners.
	 * 
	 * This method overrides the base class method to convert
	 * changes to the Job List to suitable tree changes and
	 * broadcast the changes to all registered tree listeners.
	 * This method uses the 
	 * {@link #addToPath(JobOrListWrapper, Object, ArrayList)}
	 * method to determine the tree path where the changes
	 * occurred. Then depending on the type of change
	 * it invokes the 
	 * {@link TreeModelListener#treeNodesChanged(TreeModelEvent)}, or
	 * {@link TreeModelListener#treeNodesInserted(TreeModelEvent)}, or
	 * {@link TreeModelListener#treeNodesRemoved(TreeModelEvent)}
	 * method.
	 * 
	 * @param event The workspace event that encapsulates
	 * the change that has occurred to the work space. This method
	 * only processes workspace events whose entry type is
	 * {@link EntryType#JOB_LIST}. Other changes are simply
	 * ignored by this method.
	 */
	@Override
	public void workspaceChanged(WorkspaceEvent event) {
		if (!EntryType.JOB_LIST.equals(event.getEntryType())) {
			// super.workspaceChanged(event);
			return;
		}
		// Determine the tree path that has changed starting with the root
		// and recursively build the path. The path comes in reverse order
		// from the addToPath helper method.
		ArrayList<JobOrListWrapper> pathList = new ArrayList<JobOrListWrapper>();
		addToPath(root, event.getSource(), pathList);
		// Reverse the path to obtain final TreePath.
		JobOrListWrapper path[] = new JobOrListWrapper[pathList.size()];
		for(int srcPos = pathList.size() - 1, destPos = 0; (srcPos >= 0); srcPos--, destPos++) {
			path[destPos] = pathList.get(srcPos);
		}
		// Now create an event reflecting the location of the change
		TreeModelEvent tme = new TreeModelEvent(event.getSource(), path, 
				new int[]{event.getIndexPos()}, new Object[]{event.getDeletedEntry()});
		// Dispatch event to our various tree listeners
		for(TreeModelListener tml: treeModelListeners) {
			switch (event.getOperation()) {
			case INSERT:
				tml.treeNodesInserted(tme);
				break;
			case UPDATE:
				tml.treeNodesChanged(tme);
				break;
			case DELETE:
				tml.treeNodesRemoved(tme);
			}
		}
	}

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
	 * Notify all registered listeners about removal of tree node(s).
	 * 
	 * This method can be used to broadcast notification about removal
	 * of one or more tree nodes to all registered listeners. The
	 * information about the node(s) removed are encapsulated into a
	 * suitable event for reporting to the listeners.
	 * 
	 * @param tme The tree model event that encapsulates information
	 * about the node(s) that has been removed. This parameter cannot
	 * be null and should contain suitable information. Refer to the
	 * documentation on TreeModelEvent for information about the 
	 * contents. 
	 */
	public void treeNodesRemoved(TreeModelEvent tme) {
		for(TreeModelListener tml: treeModelListeners) {
			tml.treeNodesRemoved(tme);
		}
	}

	/**
	 * Notify all registered listeners about insertion of tree node(s).
	 * 
	 * This method can be used to broadcast notification about insertion
	 * of one or more tree nodes to all registered listeners. The
	 * information about the node(s) insertion are encapsulated into a
	 * suitable event for reporting to the listeners.
	 * 
	 * @param tme The tree model event that encapsulates information
	 * about the node(s) that has been inserted. This parameter cannot
	 * be null and should contain suitable information. Refer to the
	 * documentation on TreeModelEvent for information about the 
	 * contents. 
	 */
	public void treeNodesInserted(TreeModelEvent tme) {
		for(TreeModelListener tml: treeModelListeners) {
			tml.treeNodesInserted(tme);
		}
	}

	/**
	 * A generated serialization GUID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1330539036285966008L;
}
