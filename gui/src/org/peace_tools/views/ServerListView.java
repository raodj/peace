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

import javax.swing.Icon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

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
	 */
	public ServerListView() {
		super(new BorderLayout(0, 0));
		// First create the table.
		serverTable = new JTable(new ServerListTableModel()) {
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
        toolbar = new JToolBar();
        toolbar.setFloatable(false);
        toolbar.add(Utilities.createButton("images/16x16/ServerDelete.png", 
        		null, "DeleteServer", null, 
        		"Delete the currently selected server entry", false));
        toolbar.add(Utilities.createButton("images/16x16/ServerAdd.png", 
        		null, "AddServer", null, 
        		"Add a new server entry", false));
        toolbar.add(Utilities.createButton("images/16x16/ServerInfo.png", 
        		null, "ServerInfo", null, 
        		"Obtain more information on currently selected server entry", false));
        // Add tool bar to the north
        add(toolbar, BorderLayout.NORTH);
	}
	
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
	 * An instance of the cell renderer that is used to render the
	 * first column in the table with a suitable status icon.
	 */
	private final ServerNameStatusRenderer ServerNameRenderer = new ServerNameStatusRenderer();

	/**
	 * A generated serial version ID for serialization (more
	 * realistically to keep the compiler happy). 
	 */
	private static final long serialVersionUID = 80617431851108817L;

}
