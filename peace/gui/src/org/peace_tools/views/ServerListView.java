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
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.ListSelectionModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import org.peace_tools.core.AbstractMenuHelper;
import org.peace_tools.core.MainFrame;
import org.peace_tools.core.AbstractMenuHelper.ActionType;
import org.peace_tools.data.ServerListTableModel;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

/**
 * This class provides a tabular view of the list of servers that are
 * currently configured in this work space. This table uses the 
 * ServerListTableModel class that provides the Server data from the
 *  work space in a form that is easily displayed in a table.
 */
public class ServerListView extends JPanel {
	/**
	 * The default constructor. 
	 * 
	 * The default constructor sets up the server list table and 
	 * configures the table to the default configuration. 
	 * 
	 * @param mainFrame The main frame that logically owns this job list
	 * view. The main frame is primarily used as the job listener which
	 * receives notifications on job completion. 
	 */
	public ServerListView(MainFrame mainFrame) {
		super(new BorderLayout(0, 0));
		// Save reference to main frame for future reference
		this.mainFrame = mainFrame;
		// First create the table.
		model = new ServerListTableModel();
		serverTable = new JTable(model) {
			private static final long serialVersionUID = 1430052270991586572L;
			@Override
			public TableCellRenderer getCellRenderer(int row, int column) {
				if (column == 0) {
					return ServerNameRenderer; 
				}
				return super.getCellRenderer(row, column);
			}
		};
		// Ensure table rows are not too small as our icons are 16x16
		// 19 is visual magic
		serverTable.setRowHeight(Math.max(19, serverTable.getRowHeight()));
        // Setup some column properties.
        TableColumnModel tcm = serverTable.getColumnModel();
		tcm.getColumn(0).setPreferredWidth(200);
		tcm.getColumn(1).setPreferredWidth(75);
		// Set the selection model for this table.
		serverTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        // Set some table properties
		serverTable.setBorder(null);
        serverTable.setShowHorizontalLines(true);
        serverTable.setFillsViewportHeight(true);
        serverTable.setDragEnabled(false);
        serverTable.setDropTarget(null);
        serverTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        // Place the log in a scroll pane so that servers can be scrolled
        JScrollPane scroller = new JScrollPane(serverTable);
        scroller.setBorder(null);
        scroller.getViewport().setBackground(serverTable.getBackground());
        add(scroller, BorderLayout.CENTER);
        // Create the toolbar at the top of the server list
        AbstractMenuHelper helper = mainFrame.getMenuHelper(AbstractMenuHelper.HelperType.SERVER_MENU);
        toolbar = new JToolBar();
        toolbar.setFloatable(false);
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.ADD_SERVER, false));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.CONNECTION_TEST, false));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.SHOW_MY_JOBS, false));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.SHOW_ALL_JOBS, false));
        toolbar.add(helper.getTool(AbstractMenuHelper.ActionType.REMOVE_SERVER, false));
        // Add tool bar to the north
        add(toolbar, BorderLayout.NORTH);
        // Finally create pop up menu for various job options
        createPopupMenu();
        // Add a mouse handler to trigger popups
        addMouseAdapter(serverTable);
        // Setup a table listener for the server table.
        serverTable.getSelectionModel().addListSelectionListener(helper.getListSelectionListener(serverTable));
        // Also setup the helper as a model listener as well.
		model.addTableModelListener(helper);
	}
	
	/**
	 * This is a helper method to create the pop-up menu.
	 * 
	 * This is a helper method that was introduced to streamline the
	 * code in the constructor. This method creates a popup menu
	 * that provides options for handling server entries.
	 */
	private void createPopupMenu() {
		AbstractMenuHelper helper = mainFrame.getMenuHelper(AbstractMenuHelper.HelperType.SERVER_MENU);
		popupMenu = new JPopupMenu();
		popupMenu.add(helper.getMenuItem(ActionType.ADD_SERVER, false));
		popupMenu.add(helper.getMenuItem(ActionType.CONNECTION_TEST, false));		
		popupMenu.addSeparator();
		//--------------------------------------------------
		popupMenu.add(helper.getMenuItem(ActionType.SHOW_ALL_JOBS, false));
		popupMenu.add(helper.getMenuItem(ActionType.SHOW_MY_JOBS, false));
		popupMenu.addSeparator();
		//--------------------------------------------------
		popupMenu.add(helper.getMenuItem(ActionType.REMOVE_SERVER, false));
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
	private void addMouseAdapter(JComponent list) {
		list.addMouseListener(new MouseAdapter() {
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
	 * mouse button on a item in the server table. This method checks
	 * to see if the entry is valid and if so, pops up a menu with
	 * valid operations for the selected server entry.
	 * 
	 * @param me The mouse event associated with the mouse click.
	 */
	public void handlePopup(MouseEvent me) {
		int row = serverTable.rowAtPoint(me.getPoint());
		if (model.getServer(row) != null) {
	        // Select the table entry
			serverTable.setRowSelectionInterval(row, row);
		} else {
			serverTable.clearSelection();
		}
        // Setting table entry enables or disables menu options
        // Show pop-up menu.
        popupMenu.show(serverTable, me.getX(), me.getY());
	}
	
	/**
	 * Helper method to handle double click of the mouse on a list item.
	 * 
	 * This method is invoked whenever the user double clicks on a
	 * row in the job list. This method checks to see if the
	 * entry is valid and if so, opens a tab with information about 
	 * the job running on the server.
	 * 
	 * @param me The mouse event associated with the double click.
	 */
	private void handleDoubleClick(MouseEvent me) {
		assert (me.getClickCount() == 2);
		// Obtain the row that was selected 
		int row = serverTable.rowAtPoint(me.getPoint());
		// Get the job at the given row
		Server server = model.getServer(row);
		if (server == null) {
			return;
		}
		// What do we show by default about a server?
	}
	
	/**
	 * The model that we are using to render the information in the
	 * tabular view of servers in the workspace.
	 */
	private final ServerListTableModel model;
	
	/**
	 * The toolbar that contains some commonly used tools with 
	 * the servers.
	 */
	private JToolBar toolbar;
	
	/**
	 * The actual JTable that provides a graphical view of the
	 * list of servers currently configured on this workspace.
	 */
	private JTable serverTable;
	
	/**
	 * A simple renderer that uses a JLabel to render server name and
	 * icon indicating job status.
	 * 
	 * This is an internal class that is used by the server list view
	 * to suitably display the name of the server along with a status
	 * icon. This renderer is set for the first column in the server
	 * list view.
	 */
	private class ServerNameStatusRenderer extends DefaultTableCellRenderer {
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			// Let the user class set up the renderer with defaults
			super.getTableCellRendererComponent(table, value, isSelected, hasFocus, 
					row, column);
			// Set status icon based on status.
			String statusStr  = (String) table.getValueAt(row, 1);
			statusStr         = statusStr.toUpperCase();
			// Convert status string to a number to ease icon lookup.
			int    statusCode = Server.ServerStatusType.valueOf(statusStr).ordinal();
			setIcon(ServerStatusIcons[statusCode]);
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
	private static final Icon ServerStatusIcons[] = {
		Utilities.getIcon("images/16x16/ServerInstalling.png"),
		Utilities.getIcon("images/16x16/ServerError.png"),
		Utilities.getIcon("images/16x16/ServerGood.png"),
		Utilities.getIcon("images/16x16/ServerUninstalling.png"),
		Utilities.getIcon("images/16x16/ServerError.png"),
		Utilities.getIcon("images/16x16/ServerError.png")
	};
	
	/**
	 * The main frame that logically owns this job list view.
	 */
	private final MainFrame mainFrame;
	
	/**
	 * An instance of the cell renderer that is used to render the
	 * first column in the table with a suitable status icon.
	 */
	private final ServerNameStatusRenderer ServerNameRenderer = new ServerNameStatusRenderer();
	
	/**
	 * The pop up menu that is displayed when the user left-clicks
	 * on an item in the server list view. 
	 */
	private JPopupMenu popupMenu;
	
	/**
	 * A generated serial version ID for serialization (more
	 * realistically to keep the compiler happy). 
	 */
	private static final long serialVersionUID = 80617431851108817L;
}
