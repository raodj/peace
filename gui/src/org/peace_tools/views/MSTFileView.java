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
import java.awt.Dimension;

import javax.swing.Icon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.ToolTipManager;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.TreeSelectionModel;

import org.peace_tools.core.MainFrame;
import org.peace_tools.data.EST;
import org.peace_tools.data.MSTNode;
import org.peace_tools.data.MSTTreeModel;
import org.peace_tools.generic.Utilities;

/**
 * A simple tree-based view of a Minimum Spanning Tree (MST).
 * 
 *  This is a simple view of the MST data in the form of a standard
 *  JTree.
 */
public class MSTFileView extends JPanel {
	/**
	 * The constructor creates the tree view using the data from the
	 * current work space.
	 * 
	 * @param frame
	 */
	public MSTFileView(MainFrame frame, MSTTreeModel mstModel) {
		super(new BorderLayout(0, 0));
		// Create and set up the tree model to be used.
		mstTree = new JTree((this.treeModel = mstModel));
		// Set other GUI behavioral options
		mstTree.setMinimumSize(new Dimension(50, 50));
		mstTree.setEditable(false);
		mstTree.setLargeModel(true);
		mstTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
		// Save reference to main frame.
		mainFrame = frame;
		assert(mainFrame != null);
		assert(treeModel != null);
        //Enable tool tips.
        ToolTipManager.sharedInstance().registerComponent(mstTree);
        // Setup cell renderer for better icons
        mstTree.setCellRenderer(new DataSetRenderer());
        // Add tree to the center.
        JScrollPane jsp = new JScrollPane(mstTree);
        jsp.setBorder(null);
		// Create summary information.
		JTree summaryInfo = PropertiesTreeMaker.makeProperties(treeModel.getWsEntry(), mainFrame);
		// Place the summary information and the tree-table into a split panel using
		// helper method.
		JSplitPane contentPane = 
		PropertiesTreeMaker.createPropertiesLayout("MST Information", summaryInfo, jsp, null, -1);
		// Set the content pane as the main component in this view
        add(contentPane, BorderLayout.CENTER);
	}
	
	/**
	 * The static set of icons that are repeatedly used by
	 * the custom cell renderer used by this tree view class.
	 */
    private static final Icon LeafIcons[] = {
    	Utilities.getIcon("images/16x16/ESTSmall.png"),
    };
    
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
			// Format the object value with additional information.
			MSTNode node = (MSTNode) value;
			EST est      = treeModel.getESTList().getESTs().get(node.getESTIndex());
			String fullInfo = 
				String.format("EST #%d [Metric: %d, Info: %s]", node.getESTIndex(),
					(int) node.getMetric(), est.getInfo());
        	// Let the base class do the standard processing.
            super.getTreeCellRendererComponent(tree, fullInfo, sel,
                            expanded, leaf, row, hasFocus);
            if (leaf && (value instanceof MSTNode)) {
            	setIcon(LeafIcons[0]);
            }
            return this;
        }
    }
    
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * The data model that provides the various branches and leaves
	 * of a MST to be displayed in the tree.  The tree model obtains
	 * the necessary data from the currently active work space.
	 */
	private final MSTTreeModel treeModel;
    
    /**
     * The actual JTree that provides a hierarchical view of the
     * data files currently configured on this work space.
     */
    private JTree mstTree;
	
	/**
	 * A generated serial version ID.
	 */
	private static final long serialVersionUID = 3010633830928489342L;
}
