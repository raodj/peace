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

import org.peace_tools.workspace.MSTData;

/**
 * A bridge class between a MST and a JTree.
 * 
 * This class serves as a bridge between the in-memory representation of a MST
 * (represented by a recursively defined
 * {@link org.peace_tools.data.MSTNode} set of classes). This class enables
 * reusing the data set hierarchy maintained by the
 * {@link org.peace_tools.data.MSTNode} object to display it in a JTree.
 */
public class MSTTreeModel implements TreeModel {
	/**
	 * The default constructor.
	 * 
	 * The constructor merely saves the root node of the MST for displaying
	 * it in a JTree. 
	 * 
	 * @param mst The MST data structure that is used to generate 
	 * information for this tree model.
	 * 
	 * @param estList The list of ESTs for which this MST was constructed.
	 * The EST information is used to display additional EST information
	 * in the MST display.
	 * 
	 * @param wsEntry The workspace entry corresponding to the clusters being
	 * adapted by this object.
	 */
	public MSTTreeModel(MST mst, ESTList estList, MSTData wsEntry) {
		this.mst     = mst;
		this.estList = estList;
		this.wsEntry = wsEntry;
	}
	
	/**
	 * Obtain the raw MST file data associated with this model.
	 * 
	 * @return The raw MST file data associated with this model.
	 */
	public MST getMST() { return mst; }
	
	/**
	 * Obtain the list of ESTs for which this MST was generated.
	 * 
	 * @return The list of ESTs corresponding to this MST.
	 */
	public ESTList getESTList() { return estList; }
	
	/**
	 * Adds a listener to be notified when the data associated with
	 * the data set changes.
	 * 
	 * @param tml The listener to be added to the list of listeners to
	 * be notified when the data associated with the work space changes.
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
	 * a JTree under a given node.
	 * 
	 * @param parent The object whose child is to be returned.
	 * 
	 * @param index The index of the child to be returned by this method.
	 * 
	 * @see javax.swing.tree.TreeModel#getChildCount(java.lang.Object)
	 */
	@Override
	public Object getChild(Object parent, int index) {
		if (parent instanceof MSTNode) {
			MSTNode node = (MSTNode) parent;
			return node.getChildren().get(index);
		}
		// In all other cases we don't have any children
		return null;
	}

	/**
	 * This method provides the child count for a given node.
	 * The child count is zero for leaf nodes.
	 * 
	 * @param parent The object whose child counts is to be determined. 
	 */
	@Override
	public int getChildCount(Object parent) {
		if (parent instanceof MSTNode) {
			MSTNode node = (MSTNode) parent;
			if (node.getChildren() != null) {
				return node.getChildren().size();
			}
		}
		// In all other cases we don't have any children
		return 0;
	}

	/**
	 * This method provides the zero-based logical index of a child node.
	 * The child node must have been obtained through a prior call to 
	 * the getChild() method in this class.
	 *
	 * @param parent The object whose child index is to be determined.
	 * 
	 * @param child The child node whose index is to be determined.
	 * 
	 * @return The index of the child node if found. If the child node
	 * was not found, then this method returns -1.
	 * 
	 * @see javax.swing.tree.TreeModel#getChildCount(java.lang.Object)
	 */
	@Override
	public int getIndexOfChild(Object parent, Object child) {
		if (parent instanceof MSTNode) {
			MSTNode node = (MSTNode) parent;
			return node.getChildren().indexOf(child);
		}
		// In all other cases we don't have any children
		return -1;
	}

	/**
	 * Obtain the root of the data set tree.
	 * 
	 *  @return The root node of the MST that was set when this class
	 *  was instantiated.
	 */
	@Override
	public Object getRoot() {
		return mst.getRoot();
	}

	/**
	 * Method to determine if a given entry is a leaf object. 
	 * 
	 * This method returns false, if the entry is a leaf node in the MST.
	 * In all other cases, it returns true. 
	 * 
	 * @return This method returns true if the entry is to be
	 * treated as a leaf node.
	 */
	@Override
	public boolean isLeaf(Object entry) {
		if (entry instanceof MSTNode) {
			MSTNode node = (MSTNode) entry;
			return node.isLeaf();
		}
		return true;
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
     * The only event raised by this model is TreeStructureChanged with the
     * appropriate entry to be  
     */
    protected void fireTreeStructureChanged(Object oldRoot) {
        TreeModelEvent e = new TreeModelEvent(this, new Object[] {oldRoot});
        for (TreeModelListener tml : treeModelListeners) {
            tml.treeStructureChanged(e);
        }
    }

	/**
	 * Obtain the actual workspace entry whose data is contained in this
	 * model.
	 * 
	 * A handy reference to the workspace entry from which the data for this
	 * MST tree model was actually obtained. This information can be
	 * used by "view" classes to create additional views as needed.
	 * 
	 * @return The reference to the MSTData workspace entry whose
	 * data is "modeled" by this class.
	 */
	public MSTData getWsEntry() { return wsEntry; }
	
	/**
	 * A handy reference to the workspace entry from which the data for this
	 * MST tree model was actually obtained. This information can be
	 * used by "view" classes to create additional views as needed.
	 */
	private final MSTData wsEntry;
	
	/**
	 * The MST data file that provides all the necessary information
	 * for displaying the MST.
	 * 
	 * The value is set when this class is instantiated.
	 */
	private final MST mst;
	
	/**
	 * The list of ESTs for which this MST was constructed. 
	 * The value is set when this class is instantiated. The EST
	 * information is used to display additional information about
	 * the nodes in the MST.
	 */
	private final ESTList estList;
	
    /**
     * The list of tree mode listeners that were added to this tree
     * model. This list is used to notify listeners that the tree
     * has changed.
     */
    private ArrayList<TreeModelListener> treeModelListeners =
	        new ArrayList<TreeModelListener>();
}
