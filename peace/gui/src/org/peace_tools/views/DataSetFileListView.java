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
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;

import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import org.peace_tools.core.MainFrame;
import org.peace_tools.data.DataSetTableModel;
import org.peace_tools.generic.Utilities;

/**
 * This class provides a tabular view of the list of files that are
 * currently configured in this work space. This table uses the 
 * DataSetTableModel class that provides the file data from the
 *  work space in a form that is easily displayed in a table.
 */
public class DataSetFileListView extends JPanel {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor sets up the server list table and 
	 * configures the table to the default configuration.
	 * 
	 * @param frame The main frame that ultimately owns all the views
	 */
	public DataSetFileListView(MainFrame frame) {
		super(new BorderLayout(0, 0));
		// Setup reference to main frame.
		this.mainFrame = frame;
		// First create the table.
		tableModel = new DataSetTableModel();
		fileTable = new JTable(tableModel) {
			private static final long serialVersionUID = 1430052270991586572L;
			@Override
			public TableCellRenderer getCellRenderer(int row, int column) {
				if (column == 0) {
					return FileNameRenderer; 
				}
				return super.getCellRenderer(row, column);
			}
		};
        // Setup some column properties.
        TableColumnModel tcm = fileTable.getColumnModel();
		tcm.getColumn(0).setPreferredWidth(100);
		tcm.getColumn(1).setPreferredWidth(150);
        // Set some table properties
		fileTable.setBorder(null);
        fileTable.setShowHorizontalLines(true);
        fileTable.setFillsViewportHeight(true);
        fileTable.setDragEnabled(false);
        fileTable.setDropTarget(null);
        fileTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        fileTable.setGridColor(new Color(0xd0, 0xd0, 0xd0));
        // Place the log in a scroll pane so that servers can be scrolled
        JScrollPane scroller = new JScrollPane(fileTable);
        scroller.setBorder(null);
        scroller.getViewport().setBackground(fileTable.getBackground());
        add(scroller, BorderLayout.CENTER);
        // Create the toolbar at the top of the server list
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
        // Add a double click handler 
        addMouseAdapter(fileTable);
        // Create and setup the popup menu.
        popupMenu = new DataSetPopupMenu(mainFrame, false);
	}
	
	/**
	 * A refactored helper method to add a mouse adapter. This method
	 * adds a mouse adapter to intercept certain mouse events occurring
	 * on the data set tree to trigger various operations. The mouse
	 * adapter simply delegates the actual operations to other methods 
	 * in this class.
	 * 
	 * @param list The list object to which the mouse adapter is to be
	 * added.
	 */
	private void addMouseAdapter(JComponent tree) {
		tree.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					handleDoubleClick(e);
				}
			}
			@Override
			public void mousePressed(MouseEvent e) {
				if (e.isPopupTrigger()) {
					handlePopup(e);
				}
			}
			@Override
			public void mouseReleased(MouseEvent e) {
				if (e.isPopupTrigger()) {
					handlePopup(e);
				}
			}
		});
	}
	
	/**
	 * Helper method to left mouse click on a table item.
	 * 
	 * This method is invoked whenever the user clicks on the left
	 * mouse button on a item in the table. This method checks
	 * to see if the entry is valid and if so, pops up a menu with valid
	 * operations for the selected item.
	 * 
	 * @param me The mouse event associated with the mouse click.
	 */
	public void handlePopup(MouseEvent me) {
		int row = fileTable.rowAtPoint(me.getPoint());
		// Get the file at the given row and column.
		Object entry = tableModel.getEntry(row);
		if (entry == null) {
			row = -1;
			fileTable.clearSelection();
		} else {
			// Select the table entry
			fileTable.setRowSelectionInterval(row, row);
		}
		// Update cross reference if set.
		if (treeView != null) {
			treeView.setSelectedEntry(entry);
        }
		if (entry == null) {
			return;
		}
        // Now enable/disable popup menu items based on the item
        // selected in the list.
        popupMenu.updateItems(entry, false);
        // Show pop-up menu.
        popupMenu.show(fileTable, me.getX(), me.getY());
	}
	
	/**
	 * Helper method to handle double click of the mouse on a list item.
	 * 
	 * This method is invoked whenever the user double clicks on a
	 * row in the  set tree. This method checks to see if the
	 * entry is valid and if so, opens a view of the specified
	 * data file. 
	 * 
	 * @param me The mouse event associated with the double click.
	 */
	private void handleDoubleClick(MouseEvent me) {
		assert (me.getClickCount() == 2);
		// Obtain the row that was selected 
		int row = fileTable.rowAtPoint(me.getPoint());
		// Get the file at the given row and column.
		Object entry = tableModel.getEntry(row);
		if (entry == null) {
			return;
		}
        mainFrame.getViewFactory().createView(entry, false, false);
	}
	
	/**
	 * Method to select entry in the data set table view 
	 * 
	 * This is a helper method that is primarily used by the
	 * DataSetTreeView to update the selected entry when when 
	 * the user clicks on a given entry.
	 * 
	 * @param entry The entry to be selected. If this value is null
	 * or if the entry could not be found then selections are cleared.
	 */
	public void setSelectedEntry(Object entry) {
		int row = (entry != null) ? tableModel.getRow(entry) : -1;
		if (row == fileTable.getSelectedRow()) {
			// The entry is already selected. Do nothing.
			return;
		}
		if (row == -1) {
			// No real selection.
			fileTable.clearSelection();
		} else {
			fileTable.setRowSelectionInterval(row, row);
		}
	}
	
	/**
	 * Set the reference to the data set tree view to be synchronized 
	 * with this table.
	 * 
	 * This method is used by the ViewFactory to setup cross reference
	 * between the tree view and the data set file list view.
	 * 
	 * @param treeView The tree view whose selected entry is to be 
	 * synchronized with that of this table.
	 */
	public void setDataSetTreeView(DataSetTreeView treeView) {
		this.treeView = treeView;
	}
	
	/**
	 * The tool bar that contains some commonly used tools with 
	 * the servers.
	 */
	private JToolBar toolbar;
	
	/**
	 * The actual JTable that provides a graphical view of the
	 * list of servers currently configured on this work space.
	 */
	private JTable fileTable;
	
	/**
	 * The table model that provides the data for display in the
	 * file table. This object is used to directly obtain the
	 * object at a given row.
	 */
	private DataSetTableModel tableModel;
	
	/**
	 * A inner class to serve as a special renderer that puts
	 * file names along with a suitable icon.
	 */
	private class ServerNameStatusRenderer extends JLabel
    implements TableCellRenderer {	
		/** Default constructor.
		 * 
		 * The default constructor merely sets the default font to be used
		 * by the label to be non-bold.
		 */
		public ServerNameStatusRenderer() {
			// Make the font non-bold.
			Utilities.adjustFont(this, 0, 8, -1);
			setOpaque(true);
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			setBackground(isSelected ? table.getSelectionBackground() : table.getBackground());
			// Update information.
			File file = (File) value;
			// Update the information in this JLabel
			setText(file.getName());
			if (file.exists()) {
				// Set icon based on file extension "mst", "fasta", "cls"
				String ext = Utilities.getExtension(file);
				int extType = "      fasta mst   cls".indexOf(ext) / 6;
				setIcon(FileTypeIcons[extType]);
			} else {
				// File does not exist. Set error icon.
				setIcon(FileTypeIcons[FileTypeIcons.length - 1]);
			}
			// Return the configured label for rendering.
			return this;
		}
		
		/**
		 * The general serial version UID to enable serialization of this class as
		 * per Java requirements.
		 */
		private static final long serialVersionUID = 5654108539980884223L;
	}

	/**
	 * The list of icons that are used by this table cell renderer
	 * to provide visual cues about the current status of this
	 * server.
	 */
	private static final Icon FileTypeIcons[] = {
		Utilities.getIcon("images/16x16/File.png"),
		Utilities.getIcon("images/16x16/EST.png"),
		Utilities.getIcon("images/16x16/MST.png"),
		Utilities.getIcon("images/16x16/Cluster.png"),
		Utilities.getIcon("images/16x16/Error.png")
	};
	
	/**
	 * This is a pop-up menu that is displayed whenever the user clicks
	 * the left mouse button on a table entry.  This menu is created
	 * with a static set of menu items and appropriate entries are
	 * enabled and disabled.
	 */
	DataSetPopupMenu popupMenu;

	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * An instance of the cell renderer that is used to render the
	 * first column in the table with a suitable status icon.
	 */
	private final ServerNameStatusRenderer FileNameRenderer = new ServerNameStatusRenderer();

	/**
	 * The actual tree view that provides a graphical view of the
	 * data sets in the form of a tree. This reference is used
	 * to synchronize the selected entries in the tree and 
	 * this tabular view.
	 */
    private DataSetTreeView treeView;
    
	/**
	 * A generated serial version ID for serialization (more
	 * realistically to keep the compiler happy). 
	 */
	private static final long serialVersionUID = 80617431851108817L;
}
