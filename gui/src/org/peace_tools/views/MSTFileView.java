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
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTree;
import javax.swing.ListSelectionModel;
import javax.swing.ToolTipManager;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.TreePath;
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
public class MSTFileView extends JPanel implements ListSelectionListener {
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
        
        // Now create and handle the degree information tree.
        JPanel degreeInfo = buildDegreePanel(mstModel);
        JSplitPane mstPane = 
    		PropertiesTreeMaker.createPropertiesLayout("Degree Information", degreeInfo, 
    				jsp, null, -1, null, null, null);
        
		// Create summary information.
		JTree summaryInfo = PropertiesTreeMaker.makeMSTProperties(treeModel.getWsEntry(), mainFrame);
		// Place the summary information and the tree-table into a split panel using
		// helper method.
		JSplitPane contentPane = 
		PropertiesTreeMaker.createPropertiesLayout("MST Information", summaryInfo, mstPane, null, -1);
		// Set the content pane as the main component in this view
        add(contentPane, BorderLayout.CENTER);
	}
	
	/**
	 * This method is called whenever the user selects a different row in the
	 * degree table. This method uses the EST index from the currently selected
	 * row to select the corresponding node in the adjacent tree view.
	 * 
	 * @param lse The list event associated with this call. Currently this method
	 * is ignored and the currently selected row is directly obtained from the table.
	 */
	@Override
	public void valueChanged(ListSelectionEvent lse) {
		if (degreeTable.getSelectedRow() == -1) {
			return;
		}
		// First obtain EST at the specified index
		int estIndex = (Integer) degreeTable.getValueAt(degreeTable.getSelectedRow(), 0);
		// Obtain the path to the corresponding node in the tree view.
		TreePath path = this.treeModel.getPath(estIndex);
		// Ask the tree to select the appropriate node.
		mstTree.setSelectionPath(path);
	}
	
	/**
	 * Helper method to create the table that lists the degrees of all non-leaf
	 * MST nodes. This method is called only once from the constructor. This 
	 * method was introduced to streamline the constructor and keep the code 
	 * clutter to a minimum. This method utilizes the recursive helper method
	 * {@link #populateDegreeInfo(TreeMap, MSTNode)} to obtain the degree
	 * information. It then creates a suitable JTable and populates the degree
	 * information.
	 * 
	 * @return This method returns a panel containing the degree information to
	 * be displayed. 
	 */
	private JPanel buildDegreePanel(MSTTreeModel mstModel) {
		// First get the helper method to recursively build the degree information
		// for each parent node.
		TreeMap<Integer, Integer> degreeInfo = new TreeMap<Integer, Integer>();
		populateDegreeInfo(degreeInfo, mstModel.getMST().getRoot());
		// Convert the Tree map to an 2-D array of objects to ease creation of table.
		Object[][] tableData = new Object[degreeInfo.size()][2];
		int currEntry = 0;
		for(Map.Entry<Integer, Integer> entry: degreeInfo.entrySet()) {
			tableData[currEntry][0] = entry.getKey();
			tableData[currEntry][1] = entry.getValue();
			currEntry++;
		}
		// Now create the JTable that will contain the EST index and degree info.
		final String ColTitles[] = {"EST Index", "Degree"};
		degreeTable = new JTable(tableData, ColTitles);
		// Set some of the basic properties for this degree table.
		degreeTable.setShowHorizontalLines(true);
		degreeTable.setShowVerticalLines(false);
		degreeTable.setFillsViewportHeight(true);
		degreeTable.setDragEnabled(false);
		degreeTable.setDropTarget(null);
		degreeTable.setGridColor(new Color(0xe0, 0xe0, 0xe0));
		degreeTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		degreeTable.getSelectionModel().addListSelectionListener(this);
		// Prevent reordering of columns and rows
		degreeTable.getTableHeader().setReorderingAllowed(false);
		// Setup a non-default row sorter to sort values as numbers.
		Comparator<Integer> comparator = new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return o1.compareTo(o2);
			}
		};
		// Create sorter with custom comparator
		TableRowSorter<TableModel> sorter = 
			new TableRowSorter<TableModel>(degreeTable.getModel());
		sorter.setComparator(0, comparator);
		sorter.setComparator(1, comparator);
		degreeTable.setRowSorter(sorter);
		
		// Setup column sizes
		final TableColumnModel tcm = degreeTable.getColumnModel();
		tcm.getColumn(0).setPreferredWidth(50);
		tcm.getColumn(1).setPreferredWidth(50);

		// Place the table in a scroll pane so that it can scroll
		JScrollPane scroller = new JScrollPane(degreeTable);
		scroller.getViewport().setBackground(degreeTable.getBackground());
		
		// Package the tree along with a brief formatted message into a
		// suitable panel.
		JLabel msg = new JLabel(DEGREE_TREE_INFO);
		msg.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5),
				BorderFactory.createTitledBorder("")));
		JPanel wrapper = new JPanel(new BorderLayout(5, 5));
		wrapper.add(msg, BorderLayout.NORTH);
		wrapper.add(scroller, BorderLayout.CENTER);
		// Let the constructor have the complete panel for display
		return wrapper;
	}
	
	/**
	 * This is a helper method that is used to recursively determine the degree 
	 * (number of children) for each non-leaf node in the MST.
	 * 
	 *  This is a helper method that is invoked from {@link #buildDegreePanel()} 
	 *  method to obtain the information about the degree of each non-leaf
	 *  node in the MST to be displayed in the degree table.  This method is a
	 *  recursive method that recursively descends the MST (in a depth-first)
	 *  manner and adds information about each non-leaf node to degreeInfo.
	 *   
	 * @param degreeInfo The map to be populated with the EST index &rarr; degree
	 * information. The EST index is the key into the map.
	 */
	private void populateDegreeInfo(TreeMap<Integer, Integer> degreeInfo, MSTNode node) {
		if (node.isLeaf()) {
			return;
		}
		// The node is a non-leaf node.
		final ArrayList<MSTNode> children = node.getChildren();
		// Add its degree information to the map.
		degreeInfo.put(node.getESTIndex(), children.size());
		// Now, let each of its child nodes add information about themselves
		for(MSTNode childNode: children) {
			populateDegreeInfo(degreeInfo, childNode);
		}
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
     * This table is used to display the degree (number of children) information
     * for each non-leaf node in the MST. This table is created and 
     * populated by the {@link #buildDegreePanel(MSTTreeModel)} method when
     * this class is instantiated.
     */
    private JTable degreeTable;
    
    /**
     * This instance variable contains a fixed static message that is displayed
     * to briefly explain the content and use of the degree information table
     * that is displayed in this view.
     */
    private static final String DEGREE_TREE_INFO = "<html>" +
    	"<font size=\"-1\"><i>" +
    	"The following table lists the degree of each non-leaf " +
    	"node in the MST. The entries in this table and the " +
    	"nodes in adjacent tree view are linked."    +
    	"</i></font>" +
    	"</html>";
    
	/**
	 * A generated serial version ID.
	 */
	private static final long serialVersionUID = 3010633830928489342L;
}
