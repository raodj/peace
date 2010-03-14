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

import java.util.ArrayList;

import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;

import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;
import org.peace_tools.workspace.Workspace;
import org.peace_tools.workspace.WorkspaceEvent;
import org.peace_tools.workspace.WorkspaceListener;

/**
 * A bridge class between a Workspace and a JTree.
 * 
 * This class serves as a bridge between the in-memory representation of a
 * Workspace (represented by the set of classes in the <code>workspace</code>
 * package. This class enables reusing the data set hierarchy maintained by
 * the Workspace object to display it in a JTree. In addition, this class
 * also acts to monitor and update data set views.
 */
public class DataSetTreeModel implements TreeModel, WorkspaceListener {
	/**
	 * The default constructor.
	 * 
	 * The constructor registers this object with the work space
	 * to receive notifications on Job status updates.
	 */
	public DataSetTreeModel() {
		Workspace.get().addWorkspaceListener(this);
	}
	
	/**
	 * Adds a listener to be notified when the data associated with
	 * the data set changes.
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
	 * This method provides the actual child object to be displayed in
	 * a JTree depending on the data in the Workspace. Specifically this
	 * method returns the following object depending on the data type 
	 * of parent:
	 * 
	 *  <ul>
	 *  
	 *  <li>if parent is the globally unique workspace instance, then
	 *  this method returns the DataSet object at the given index.</li>
	 *  
	 *  <li>if parent is a DataSet, then this method returns either an
	 *  MSTData object (if index < number of MST entries in DataSet) or
	 *  a MSTClusterData object.</li>
	 *  
	 *  <li>In all other cases, this method returns null.</li>
	 *  
	 *  </ul>
	 * 
	 * @param parent The object whose child is to be returned.
	 * 
	 * @param index The index of the child to be returned by this method.
	 * 
	 * @see javax.swing.tree.TreeModel#getChildCount(java.lang.Object)
	 */
	@Override
	public Object getChild(Object parent, int index) {
		if (parent instanceof Workspace) {
			assert ( parent == Workspace.get() );
			Workspace root = (Workspace) parent;
			return root.getDataSets().get(index);
		} else if (parent instanceof DataSet) {
			DataSet ds = (DataSet) parent;
			if (index < ds.getMSTList().size()) {
				return ds.getMSTList().get(index);
			}
			index -= ds.getMSTList().size();
			return ds.getClusterList().get(index);
		}
		// In all other cases we don't have any children
		return null;
	}

	/**
	 * This method provides the child count depending on the data 
	 * in the Workspace. Specifically this method returns the following
	 * child counts depending on the data type of parent:
	 * 
	 *  <ul>
	 *  
	 *  <li>if parent is the globally unique workspace instance, then
	 *  this method returns the number of data sets in the work space
	 *  as the number of children.</li>
	 *  
	 *  <li>if parent is a DataSet, then this method returns the number
	 *  of MST and Cluster data files contained in the data set.</li>
	 *  
	 *  <li>In all other cases, this method returns zero.</li>
	 *  
	 *  </ul>
	 * 
	 * @param parent The object whose child counts is to be determined. 
	 * 
	 * @see javax.swing.tree.TreeModel#getChildCount(java.lang.Object)
	 */
	@Override
	public int getChildCount(Object parent) {
		if (parent instanceof Workspace) {
			assert ( parent == Workspace.get() );
			Workspace root = (Workspace) parent;
			return root.getDataSets().size();
		} else if (parent instanceof DataSet) {
			DataSet ds = (DataSet) parent;
			return ds.getMSTList().size() + ds.getClusterList().size();
		}
		// In all other cases we don't have any children
		return 0;
	}

	/**
	 * This method provides the zero-based logical index of a child node.
	 * The child node must have been obtained through a prior call to 
	 * the getChild() method in this class.
	 * 
	 * Specifically this method returns the following
	 * child counts depending on the data type of parent:
	 * 
	 *  <ul>
	 *  
	 *  <li>if parent is the globally unique workspace instance, then
	 *  this method assumes child is a DataSet and returns the index of
	 *  the child in the list of data sets in the workspace.</li>
	 *  
	 *  <li>if parent is a DataSet, then this method assumes that the child
	 *  is either a MSTData or MSTClusterData object and returns index from
	 *  either of the two lists.</li>returns the number
	 *  
	 *  <li>In all other cases, this method returns -1.</li>
	 *  
	 *  </ul>
	 * 
	 * @param parent The object whose child counts is to be determined. 
	 * 
	 * @see javax.swing.tree.TreeModel#getChildCount(java.lang.Object)
	 */
	@Override
	public int getIndexOfChild(Object parent, Object child) {
		if (parent instanceof Workspace) {
			assert ( parent == Workspace.get() );
			Workspace root = (Workspace) parent;
			assert(child instanceof DataSet);
			return root.getDataSets().indexOf(child);
		} else if (parent instanceof DataSet) {
			DataSet ds = (DataSet) parent;
			assert((child instanceof MSTData) || (child instanceof MSTClusterData));
			int index = ds.getMSTList().indexOf(child);
			if (index == -1) {
				index = ds.getClusterList().indexOf(child);
				if (index > -1) {
					index += ds.getMSTList().size();
				}
			}
			return index;
		}
		// In all other cases we don't have any children
		return -1;
	}

	/**
	 * Obtain the root of the data set tree. The workspace entry is 
	 * simply used as the root of the JTree tree.
	 */
	@Override
	public Object getRoot() {
		Workspace workspace = Workspace.get();
		return workspace;
	}

	/**
	 * Method to determine if a given entry is a leaf object. 
	 * 
	 * This method returns false, if the entry is a Workspace or
	 * a DataSet. In all other cases, it returns true. 
	 * 
	 * @return This method returns true if the entry is to be
	 * treated as a leaf node.
	 */
	@Override
	public boolean isLeaf(Object entry) {
		return !((entry instanceof Workspace) || 
				(entry instanceof DataSet));
	}

	/**
	 * This data model is meant to be used with a read-only type
	 * JTree. Consequently, this method (that is meant to change
	 * the value in the tree) is not implemented and must not
	 * be called.
	 */
	@Override
	public void valueForPathChanged(TreePath arg0, Object arg1) {
		assert("Not implemented" == null);
	}

	/**
	 * This is a helper method to determine path to a given entry.
	 * 
	 * This method can be used to obtain the path within the 
	 * data set for a given data set entry. This is 
	 *  
	 * @param entry The entry for which the tree path is required.
	 * 
	 * @return The complete path to the specified entry. If the
	 * specified path could not be determined then this method
	 * returns null.
	 */
	public TreePath getPath(Object entry) {
		if (entry instanceof Workspace) {
			// This is the root entry. Simple case.
			return new TreePath(entry);
		}
		if ((entry instanceof DataSet) && 
			(getIndexOfChild(getRoot(), entry) != -1)) {
			// This data set is a valid entry.
			return new TreePath(new Object[]{getRoot(), entry});
		}
		if (entry instanceof MSTData) {
			// Return the path to the MST entry
			MSTData mst = (MSTData) entry;
			return new TreePath(new Object[]{getRoot(), mst.getDataSet(), mst});
		}
		if (entry instanceof MSTClusterData) {
			// Return the path to the MST entry
			MSTClusterData cls = (MSTClusterData) entry;
			return new TreePath(new Object[]{getRoot(), cls.getDataSet(), cls});
		}
		// Can't figure out the path
		return null;
	}
	
    /**
     * The only event raised by this model is TreeStructureChanged with the
     * appropriate entry to be  
     */
    protected void fireTreeStructureChanged(Object oldRoot) {
        TreeModelEvent e = new TreeModelEvent(this, new Object[] {oldRoot});
        for (TreeModelListener tml : treeModelListeners) {
            tml.treeStructureChanged(e);
        }
    }

	@Override
	public void workspaceChanged(WorkspaceEvent event) {
		// This data set tree model only cares about events relating to data set, mst
		// and cluster information. Everything else is ignored.
		final WorkspaceEvent.EntryType entryType = event.getEntryType();
		if (!entryType.equals(WorkspaceEvent.EntryType.DATA_SET) &&
			!entryType.equals(WorkspaceEvent.EntryType.MST_DATA) &&
			!entryType.equals(WorkspaceEvent.EntryType.MST_CLUSTER_DATA)) {
			// We don't care bout this update.
			return;
		}
		// Translate work space event to a model event.
		if (event.getOperation().equals(WorkspaceEvent.Operation.UPDATE)) {
			fireTreeStructureChanged(event.getSource());
			return;
		} 
		// For inserts and deletes the changes are a bit more
		// extensive to the tree.
		if (event.getEntryType().equals(WorkspaceEvent.EntryType.DATA_SET)) {
			// This event is at the data set level.
			fireTreeStructureChanged(Workspace.get());
		} else {
			// This event is within a data set. Fire change at that level
			Object oldRoot = null;
			if (event.getEntryType().equals(WorkspaceEvent.EntryType.MST_DATA)) {
				oldRoot = ((MSTData) event.getSource()).getDataSet();
			} else if (event.getEntryType().equals(WorkspaceEvent.EntryType.MST_CLUSTER_DATA)) {
				oldRoot = ((MSTClusterData) event.getSource()).getDataSet();
			}
			if (oldRoot != null) {
				fireTreeStructureChanged(oldRoot);
			}
		}
	}
	
    /**
     * The list of tree mode listeners that were added to this tree
     * model. This list is used to notify listeners that the tree
     * has changed.
     */
    private ArrayList<TreeModelListener> treeModelListeners =
	        new ArrayList<TreeModelListener>();
}
