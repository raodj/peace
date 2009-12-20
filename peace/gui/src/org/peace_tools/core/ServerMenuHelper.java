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

import javax.swing.Box;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JToolBar;

import org.peace_tools.core.server.ServerWizard;
import org.peace_tools.generic.Utilities;

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
public class ServerMenuHelper implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public ServerMenuHelper(MainFrame mainFrame) {
		this.mainFrame = mainFrame;
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
		JMenuItem item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Add New Server",
					"Launch wizard to add a new server and install PEACE runtime on server.",
					"AddServer", this, "images/24x24/ServerAdd.png", 
					null, true, false);
		serverMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Remove Server Entry",
					"Uninstall PEACE runtime and remove an existing server entry from he workspace.", 
					"RemoveServer", this, "images/24x24/ServerDelete.png", 
					null, true, false);
		serverMenu.add(item);
		
		serverMenu.addSeparator();
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"View All Servers",
				"A tabular view of all the servers currently configured in this workspace",
				"ViewServers", this, "images/24x24/Server.png", null, true, false);
		serverMenu.add(item);
		
		// Add tool bar entry to add a server entry.
		if (toolbar != null) {
			// Add some of the tools that we anticipate users to work
			// with frequently to the tool bar.
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(Utilities.createToolButton("images/24x24/ServerAdd.png", 
					null, "AddServer", this,
					"Launch wizard to add a new server and install " +
					"PEACE runtime on server.", true));
		}
		
		return serverMenu;
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
		String command = event.getActionCommand();
		if ("AddServer".equals(command)) {
			addServer();
		}
	}

	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
}
