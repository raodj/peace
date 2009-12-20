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

import org.peace_tools.generic.Utilities;

/**
 * The help menu helper for the application.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "View" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that the ViewMenuHelper is essentially an event
 * handler that is set on the various menu items in this class. This helper
 * class also provides a  createViewMenu() method that actually creates
 * the "View" menu.
 */
public class ViewMenuHelper implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public ViewMenuHelper(MainFrame mainFrame) {
		this.mainFrame = mainFrame;
		assert(this.mainFrame != null);
	}

	/**
	 * Helper method to create the view menu.
	 * 
	 * This method performs the actual task of creating the "View" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "View" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "View" menu.
	 * 
	 * @param toolbar The tool bar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the tool bar is null, then
	 * shortcuts are not added.
	 */
	public JMenu createViewMenu(JToolBar toolbar) {
		// Create the actual menu.
		JMenu viewMenu = new JMenu("View  ");
	    // First create and add the new file creation menu options.
		JMenuItem item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "New Main Frame",
				"Open a new main frame to provide extra screen real estate to view data",
				"NewFrame", this, "images/24x24/MainFrame.png", null, true, false);
		viewMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, "Duplicate Editor",
				"Add a duplicate view of the currently selected view to the same main frame",
				"DupTab", this, "images/24x24/DupTab.png", 
				null, true, false);
		viewMenu.add(item);
		viewMenu.addSeparator();
		//----------------------------------------------------------
		// Menu items to display any special views
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Workspace Hierarchy Browser",
				"View files in workspace in a hierarchy based on data sets",
				"DATASET_TREE", this, "images/24x24/DataSetView.png", null, true, false);
		viewMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Workspace File Browser",
				"View files in current workspace in as a simple but comprehensive list",
				"DATASET_FILE", this, "images/24x24/FileView.png", null, true, false);
		viewMenu.add(item);
		viewMenu.addSeparator();
		//----------------------------------------------------------
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"View All Jobs",
				"View the current and completed jobs in the current workspace",
				"JOB_LIST", this, "images/24x24/Jobs.png", null, true, false);
		viewMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"View All Servers",
				"A tabular view of all the servers currently configured in this workspace",
				"SERVER_LIST", this, "images/24x24/Server.png", null, true, false);
		viewMenu.add(item);
		viewMenu.addSeparator();
		//----------------------------------------------------------
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"User Logs",
				"Opens a tabular view of logs that provide additional information to a user",
				"USER_LOGS", this, "images/24x24/UserLog.png", 
				null, true, false);
		viewMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Programmer Logs",
				"View a textual form of the programmer logs with internal logs",
				"PROGRAMMER_LOGS", this, "images/24x24/ProgLog.png", 
				null, true, false);
		viewMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"Welcome & Quick Start",
				"The initial welcome screen with quick start guide about PEACE",
				"WELCOME_SCREEN", this, "images/24x24/WelcomeSrc.png", 
				null, true, false);
		viewMenu.add(item);
		
		if (toolbar != null) {
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(Utilities.createToolButton("images/24x24/UserLog.png", 
					null, "ViewUserLogs", this,
					"Opens a tabular view of logs that provide additional information to a user", 
					true));
		}
		return viewMenu;
	}
	
 
	@Override
	public void actionPerformed(ActionEvent event) {
		String cmd = event.getActionCommand();
		// Process the command 
		if ("NewFrame".equals(cmd)) {
			// Not yet implemented
		} else if ("DupTab".equals(cmd)) {
			// Not yet implement
		} else if ("WELCOME_SCREEN".equals(cmd)) {
			// Show the welcome screen if it is not already there.
			ViewFactory vf = mainFrame.getViewFactory();
			vf.createView("installFiles/welcome.html", null,
						  ViewFactory.ViewType.HTML_VIEW, false, false);
		} else {
			// This must be to view a specific view. The action
			// command has the type of view to be shown.
			ViewFactory.ViewType view = ViewFactory.ViewType.valueOf(cmd);
			if (view != null) {
				// Display the static view if it is not already visible.
				mainFrame.getViewFactory().showStaticView(view);
			}
		}
	}

	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
}
