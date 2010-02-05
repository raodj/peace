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

package org.peace_tools.views;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.Icon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.ToolTipManager;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.peace_tools.core.MainFrame;
import org.peace_tools.data.DataSetTreeModel;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;

public class DataSetTreeView extends JPanel implements ActionListener {
	/**
	 * The constructor creates the tree view using the data from the
	 * current work space.
	 * 
	 * @param frame The main frame that logically owns this view.
	 */
	public DataSetTreeView(MainFrame frame) {
		super(new BorderLayout(0, 0));
		// Create and set up the tree model to be used.
		treeModel   = new DataSetTreeModel();
		dataSetTree = new JTree(treeModel);
		// Set other GUI behavioral options
		dataSetTree.setEditable(false);
		dataSetTree.setLargeModel(true);
		dataSetTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
		// Ensure each row in tree is slightly bigger than our icon.
		dataSetTree.setRowHeight(Math.max(19, dataSetTree.getRowHeight()));
		// Save reference to main frame.
		mainFrame = frame;
		assert(mainFrame != null);
        //Enable tool tips.
        ToolTipManager.sharedInstance().registerComponent(dataSetTree);
        // Setup cell renderer for better icons
        dataSetTree.setCellRenderer(new DataSetRenderer());
        // Add a mouse adapter to intercept double clicks to open
        // the files on the tree entry items.
        addMouseAdapter(dataSetTree);
        // Wrap the tree view in a scroll pane to let it scroll.
		JScrollPane jsp = new JScrollPane(dataSetTree);
		jsp.setBorder(null);
        this.add(jsp, BorderLayout.CENTER);
        
        // Create and setup the tool bar.
        toolbar = new JToolBar();
        toolbar.setFloatable(false);
        toolbar.add(Utilities.createButton("images/16x16/MST.png", 
                        null, "NewMST", null, 
                        "Compute MST & Clusters for selected data set", false));
        toolbar.add(Utilities.createButton("images/16x16/Cluster.png", 
                        null, "NewClusters", null, 
                        "Compute clusters based on selected MST", false));
        // Add tool bar to the north
        add(toolbar, BorderLayout.NORTH);
        
        // Create popup menu for later display
        popupMenu = new DataSetPopupMenu(mainFrame, true);
        popupMenu.addActionListener(this);
	}

	/**
	 * A refactored helper method to add a mouse adapter. This method
	 * adds a mouse adapter to intercept certain mouse events occurring
	 * on the data set tree to trigger various operations. The mouse
	 * adapter simply delegates the actual operations to other methods 
	 * in this class.
	 * 
	 * @param tree The tree object to which the mouse adapter is to be
	 * added.
	 */
	private void addMouseAdapter(JTree tree) {
		tree.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					handleDoubleClick(e);
				}
			}
			@Override
			public void mousePressed(MouseEvent e) {
				// Needed for popup menus under Linux
				if (e.isPopupTrigger()) {
					handlePopup(e);
				}
			}
			@Override
			public void mouseReleased(MouseEvent e) {
				// Needed for popup menus under Windows
				if (e.isPopupTrigger()) {
					handlePopup(e);
				}
			}
		});
	}
	
	/**
	 * Helper method to left mouse click on a tree item.
	 * 
	 * This method is invoked whenever the user clicks on the left
	 * mouse button on a item in the data set tree. This method checks
	 * to see if the entry is valid and if so, pops up a menu with valid
	 * operations for the selected item.
	 * 
	 * @param me The mouse event associated with the mouse click.
	 */
	public void handlePopup(MouseEvent me) {
        TreePath path = dataSetTree.getPathForLocation(me.getX(), me.getY());
        if (path == null) {
        	// There were no nodes at the click point. Nothing to do
        	return;
        }
        // Select the path first
        dataSetTree.setSelectionPath(path);
        // Update cross reference if set.
        if (listView != null) {
        	listView.setSelectedEntry(path.getLastPathComponent());
        }
        // Obtain needed information for pop up menu
        final Object  entry      = path.getLastPathComponent();
        final boolean isExpanded = dataSetTree.isExpanded(path);
        // Now enable/disable popup menu items based on the item
        // selected in the tree.
        popupMenu.updateItems(entry, isExpanded);
        // Show pop-up menu.
        popupMenu.show(dataSetTree, me.getX(), me.getY());
	}
	
	/**
	 * Helper method to handle double click of the mouse on a tree item.
	 * 
	 * This method is invoked whenever the user double clicks on a
	 * item in the data set tree. This method checks to see if the
	 * entry is valid and if so, opens a view of the specified
	 * data file. 
	 * 
	 * @param me The mouse event associated with the double click.
	 */
	private void handleDoubleClick(MouseEvent me) {
		assert (me.getClickCount() == 2);
		// First get the node corresponding to the location of mouse
        TreePath path = dataSetTree.getSelectionPath();
        if ((path == null) || (path.getPathCount() < 2)){
        	// There were no nodes at the click point. Nothing to do
        	return;
        }
        // Check and display the view (if one does not exist)
        mainFrame.getViewFactory().createView(path.getLastPathComponent(), false, false);
	}
	    
	/**
	 * A tree cell renderer that essentially provides a better 
	 * representative set of Icons to make the overall display 
	 * look a bit prettier.
	 */
    private class DataSetRenderer extends DefaultTreeCellRenderer {
    	private static final long serialVersionUID = -5720888296480938746L;
		@Override
        public Component getTreeCellRendererComponent(JTree tree,
        		Object value, boolean sel, boolean expanded,
                boolean leaf, int row, boolean hasFocus) {
        	// Let the base class do the standard processing.
            super.getTreeCellRendererComponent(tree, value, sel,
                            expanded, leaf, row, hasFocus);
            if (!leaf && (value instanceof DataSet)) {
            	setIcon(LeafIcons[2]);
            	setToolTipText(((DataSet) value).getDescription());
            }
            // Replace icons as needed.
            if (leaf && (value instanceof MSTData)) {
            	MSTData mst = (MSTData) value;
                setIcon(LeafIcons[0]);
                setToolTipText(mst.getDescription());
            }
            if (leaf && (value instanceof MSTClusterData)) {
            	MSTClusterData cluster = (MSTClusterData) value;
                setIcon(LeafIcons[1]);
                setToolTipText(cluster.getDescription());
            }
            return this;
        }
    }    

	/**
	 * Action listener to handle call backs from pop up menu items.
	 * 
	 * This action handler is invoked whenever the user selects various
	 * menu options or clicks on the tool bar buttons. This method is
	 * handles only actions that are specific to tree operations
	 * (such as expand and collapse). Other operations are handled
	 * by the DataSetPopupMenu class. 
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	
	/**
	 * Set the reference to the file table to be synchronized with
	 * this tree.
	 * 
	 * This method is used by the ViewFactory to setup cross reference
	 * between the tree view and the data set file list view.
	 * 
	 * @param listView The list view whose selected entry is to be 
	 * synchronized with that of this tree.
	 */
	public void setTableView(DataSetFileListView listView) {
		this.listView = listView;
	}
	
	/**
	 * Method to select entry in the data set tree view 
	 * 
	 * This is a helper method that is primarily used by the
	 * DataSetFileListView to update the selected entry when when
	 * the selection in the table changes.
	 * 
	 * @param entry The entry to be selected. If this value is null
	 * then selections are cleared.
	 */
	public void setSelectedEntry(Object entry) {
		TreePath path = (entry != null) ? treeModel.getPath(entry) : null;
		if (path == null) {
			// No real selection.
			dataSetTree.clearSelection();
		} else {
			dataSetTree.setSelectionPath(path);
			dataSetTree.expandPath(path);
		}
	}
	
	/**
	 * The static set of icons that are repeatedly used by
	 * the custom cell renderer used by this tree view class.
	 */
    private static final Icon LeafIcons[] = {
    	Utilities.getIcon("images/16x16/MST.png"),
    	Utilities.getIcon("images/16x16/Cluster.png"),
    	Utilities.getIcon("images/16x16/DataSet.png")
    };

	/**
	 * This is a pop-up menu that is displayed whenever the user clicks
	 * the left mouse button on a menu item.  This menu is created
	 * with a static set of menu items and appropriate entries are
	 * enabled and disabled.
	 */
	DataSetPopupMenu popupMenu;
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this view in its frame. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * The data model that provides the various branches and leaves
	 * to be displayed in the tree.  The tree model obtains the 
	 * necessary data from the currently active workspace.
	 */
	private final DataSetTreeModel treeModel;
	
    /**
     * The tool bar that contains some commonly used tools with 
     * the data set tree view.
     */
    private JToolBar toolbar;
    
    /**
     * The actual JTree that provides a hierarchical view of the
     * data files currently configured on this work space.
     */
    private JTree dataSetTree;
	
	/**
	 * The actual JTable that provides a graphical view of the
	 * data sets in the form of a list. This reference is used
	 * to synchronize the selected entries in the tree and 
	 * the table.
	 */
    private DataSetFileListView listView;
    
	/**
	 * A generated serial version ID.
	 */
	private static final long serialVersionUID = 3010633830928489342L;
}
