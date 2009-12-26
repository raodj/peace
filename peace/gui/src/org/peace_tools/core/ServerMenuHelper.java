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

package org.peace_tools.core;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.AbstractButton;
import javax.swing.Box;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionListener;

import org.peace_tools.core.server.ServerWizard;
import org.peace_tools.data.ServerListTableModel;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

/**
 * The server menu helper for the application.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "Server" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that the ServerMenuHelper is essentially an event
 * handler that is set on the various menu items in this class. This helper
 * class also provides a  createServerMenu() method that actually creates
 * the "Server" menu.
 */
public class ServerMenuHelper extends AbstractMenuHelper
implements ActionListener, ListSelectionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public ServerMenuHelper(MainFrame mainFrame) {
		super(HelperType.SERVER_MENU, mainFrame);
	}

	/**
	 * Helper method to create the server menu.
	 * 
	 * This method performs the actual task of creating the "Server" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "Server" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "Server" menu.
	 * 
	 * @param toolbar The toolbar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the toolbar is null, then
	 * shortcuts are not added.
	 */
	public JMenu createServerMenu(JToolBar toolbar) {
		// Create the actual menu.
		JMenu serverMenu = new JMenu("Server  ");
	    // First create and add the new file creation menu options.
		serverMenu.add(getMenuItem(ActionType.ADD_SERVER, true));
		serverMenu.addSeparator();
		serverMenu.add(getMenuItem(ActionType.SHOW_MY_JOBS, true));
		serverMenu.add(getMenuItem(ActionType.SHOW_ALL_JOBS, true));
		serverMenu.add(getMenuItem(ActionType.CONNECTION_TEST, true));
		
		serverMenu.addSeparator();
		serverMenu.add(getMenuItem(ActionType.REMOVE_SERVER, true));
		
		// Add tool bar entry to add a server entry.
		if (toolbar != null) {
			// Add some of the tools that we anticipate users to work
			// with frequently to the tool bar.
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(getTool(ActionType.ADD_SERVER, true));
		}
		
		return serverMenu;
	}
	
	@Override
	public ActionListener getActionListener() {
		return this;
	}

	@Override
	public ListSelectionListener getListSelectionListener(JTable table) {
		this.table = table;
		return this;
	}

    /**
     * The selection listener/handler for a table.
     * 
     * This method is invoked by the core Swing classes whenever the
     * user selects a specific entry in a list or table. This
     * method essentially enables/disables various tool bar buttons
     * and menu items based on the current job selection.
     * 
     * @note Currently we only handle JTable and not JList.
     * 
     * @param event The selection event associated with this method.
     * This event is not really used.
     */
	@Override
	public void valueChanged(ListSelectionEvent event) {
		if (event.getValueIsAdjusting()) {
			return;
		}
		assert ( table != null );
		if (!(table.getModel() instanceof ServerListTableModel)) {
			// The table model is not really compatible.
			return;
		}
		ServerListTableModel model = (ServerListTableModel) table.getModel();
		final int row           = table.getSelectedRow();
		// A couple of flags to help make checks below easier
		server = model.getServer(row);

		setEnabled("RemoveServer",  server != null); 
		setEnabled("ShowMyJobs",    server != null);
		setEnabled("ServerInfo",    server != null);
		setEnabled("ServerConnect", server != null);
	}
	
	@Override
	public JMenuItem getMenuItem(ActionType actionType, boolean mainMenu) {
		int index = actionType.ordinal() - ActionType.ADD_SERVER.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainMenu ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";
		// Create and return the main menu item
		JMenuItem item =
			Utilities.createMenuItem(Utilities.MENU_ITEM, MenuTitles[index],
				(mainMenu ? MenuSubTitles[index] : null),
				ActionCmds[index], this, IconPath, 
				null, true, false);
		if (index > 0) {
			// Track context sensitive entries
			contextItemList.add(item);
			Utilities.setEnabled(item, false);
		}
		return item;
	}
	
	@Override
	public AbstractButton getTool(ActionType actionType, boolean mainToolBar) {
		int index = actionType.ordinal() - ActionType.ADD_SERVER.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainToolBar ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";

		AbstractButton item =
			Utilities.createToolButton(IconPath, null, ActionCmds[index], this, 
				MenuSubTitles[index], true);
		if (index > 0) {
			// Track context sensitive entries
			contextItemList.add(item);
			item.setEnabled(false);
		}
		return item;
	}
	
    /** Add a new server entry.
     * 
     * This method is invoked when the user chooses to add a new
     * server entry. This method is actually invoked 
     * from the actionPerformed() method.  This method starts up
     * a new instance of the ServerWizard that guides the user
     * through the process of adding a new server entry and 
     * installing PEACE on the server.
     */
	 private void addServer() {
		 ServerWizard sw = new ServerWizard("Add new server", mainFrame,
				 mainFrame.getCenterPane());
		 sw.showWizard(null);
	 }

	@Override
	public void actionPerformed(ActionEvent event) {
		String cmd = event.getActionCommand();
		if ("AddServer".equals(cmd)) {
			addServer();
		} else if ("RemoveServer".equals(cmd)) {
			// Use helper method to verify and delete entry.
			deleteServer();
		} else if ("ShowMyJobs".equals(cmd)) {
			// List just the user's jobs on the server.
			mainFrame.getViewFactory().createView(server, server.getUserID());
		} else if ("ServerInfo".equals(cmd)) {
			// List all jobs on the server.
			mainFrame.getViewFactory().createView(server, null);
		} else if ("ServerConnect".equals(cmd)) {
			// Test connection.
		}
	}

	/**
	 * Helper method to setup dialog to display jobs/processes on
	 * a given server.
	 * 
	 * This method is a helper method that is used by both JobMenuHelper
	 * and ServerMenuHelper to display jobs/processes running on a
	 * given server.
	 * 
	 * @param server The server whose jobs are to be listed.
	 * @param allJobs If this flag is true then all the jobs running/
	 * scheduled on the server are displayed rather than just 
	 */
	public void showJobs(Server server, boolean allJobs) {
		
	}
	
	/**
	 * Helper method to verify and delete server entry.
	 * 
	 * This is a helper method that was introduced to keep the code
	 * clutter in the actionPerformed() method to a minimum. This
	 * method uses the DeleteDialog helper dialog to actually 
	 * delete the server entry.
	 * 
	 * @note Deleting the server entry from the work space will
	 * cause the GUI to automatically reflect the changes.
	 */
	private void deleteServer() {
		if (server == null) {
			return;
		}
		// Delete the entry via the delete dialog.
		DeleteDialog delDialog = new DeleteDialog(mainFrame, server);
		delDialog.setVisible(true);
	}
	
	@Override
	public TreeSelectionListener getTreeSelectionListener(JTree tree) {
		return null;
	}
	
	/**
	 * The currently selected server entry (if any). This value is
	 * updated by the list selection listener whenever a valid
	 * server entry is selected. If a valid server entry is not selected
	 * then this entry is set to null.
	 */
	private Server server;
	
	/**
	 * The strings for each menu title created by this helper. The list of
	 * values are organized in the same order as the ordinal values of 
	 * the ActionType enumeration. Preserving the order is important.
	 */
	private static final String MenuTitles[] = {
		"Add New Server", 
		"Remove Server Entry",
		"Show your jobs on server",
		"Show all jobs on server",
		"Test server connection"
	};
	
	/**
	 * The icon file names for each menu title created by this helper. 
	 * The list of values are organized in the same order as the ordinal 
	 * values of the ActionType enumeration. Preserving the order is
	 * important. Note that a prefix directory (such as: images/16x16)
	 * and a suffix extension (.png) is added when tools or menu items
	 * are created.
	 */
	private static final String IconNames[] = {
		"ServerAdd", 
		"ServerDelete",
		"ServerMyJobs",
		"ServerInfo",
		"ServerConnect"
	};
	
	/**
	 * The strings for the action commands generated by the various menu items
	 * in the help menu. The list of values are organized in the same order as
	 * the ordinal values of the ActionType enumeration. Preserving the order 
	 * is important.
	 */
	private static final String ActionCmds[] = {
		"AddServer", 
		"RemoveServer",
		"ShowMyJobs",
		"ServerInfo",
		"ServerConnect"
	};
	
	/**
	 * The various sub menu titles that are used in the main menu. The
	 * sub menu titles are used to provide the user with a bit more 
	 * verbose description on the action that will be performed by a given
	 * menu item. The list of values are organized in the same order as 
	 *  the ordinal values of the ActionType enumeration. Preserving the 
	 *  order is important.
	 */
	private static final String MenuSubTitles[] = {
		"Launch wizard to add a new server and install PEACE runtime on server.", 
		"Uninstall PEACE runtime and remove an existing server entry from he workspace.",
		"Show just my jobs that are running or queued on the server",
		"Show all the jobs that are running or queued on the server",
		"Try connecting to the server to ensure communication is operational"
	};
}
