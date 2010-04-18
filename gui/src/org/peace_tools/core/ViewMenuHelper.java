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
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionListener;

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
public class ViewMenuHelper extends AbstractMenuHelper 
implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "Help" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public ViewMenuHelper(MainFrame mainFrame) {
		super(HelperType.VIEW_MENU, mainFrame);
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
	 * 
	 * @param hmh The help menu helper object associated with the main 
	 * frame. It is used to add a menu entry for the quick start guide.
	 */
	public JMenu createViewMenu(JToolBar toolbar, AbstractMenuHelper hmh) {
		// Create the actual menu.
		JMenu viewMenu = new JMenu("View  ");
	    // First create and add the new file creation menu options.
		viewMenu.add(getMenuItem(ActionType.NEW_MAIN_FRAME, true));
		viewMenu.add(getMenuItem(ActionType.DUPLICATE_EDITOR, true));
		viewMenu.addSeparator();
		//----------------------------------------------------------
		// Menu items to display any special views
		viewMenu.add(getMenuItem(ActionType.SHOW_WORKSPACE_HIERARCHY_BROWSER, true));
		viewMenu.add(getMenuItem(ActionType.SHOW_WORKSPACE_FILE_BROWSER, true));
		viewMenu.addSeparator();
		//----------------------------------------------------------
		viewMenu.add(getMenuItem(ActionType.VIEW_JOBS, true));
		viewMenu.add(getMenuItem(ActionType.VIEW_SERVERS, true));
		viewMenu.addSeparator();
		//----------------------------------------------------------
		viewMenu.add(getMenuItem(ActionType.SHOW_USER_LOGS, true));
		viewMenu.add(getMenuItem(ActionType.SHOW_PROG_LOGS, true));
		viewMenu.add(hmh.getMenuItem(ActionType.SHOW_QUICK_START, true));
		
		if (toolbar != null) {
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(getTool(ActionType.SHOW_USER_LOGS, true));
		}
		return viewMenu;
	}
	
	@Override
	public JMenuItem getMenuItem(ActionType actionType, boolean mainMenu) {
		int index = actionType.ordinal() - ActionType.NEW_MAIN_FRAME.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainMenu ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";
		// Create and return the main menu item
		return Utilities.createMenuItem(Utilities.MENU_ITEM, MenuTitles[index],
				(mainMenu ? MenuSubTitles[index] : null),
				ActionCmds[index], this, IconPath, 
				null, true, false);
	}
	
	@Override
	public AbstractButton getTool(ActionType actionType, boolean mainToolBar) {
		int index = actionType.ordinal() - ActionType.NEW_MAIN_FRAME.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainToolBar ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";

		return Utilities.createToolButton(IconPath, null, ActionCmds[index], this, 
				MenuSubTitles[index], true);
	}
	
	@Override
	public TreeSelectionListener getTreeSelectionListener(JTree tree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ActionListener getActionListener() {
		return this;
	}

	@Override
	public ListSelectionListener getListSelectionListener(JTable table) {
		return null;
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
			vf.createView("installFiles/welcome.html", null, null,
						  ViewFactory.ViewType.HTML_VIEW, false, false, null);
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
	 * The strings for each menu title created by this helper. The list of
	 * values are organized in the same order as the ordinal values of 
	 * the ActionType enumeration. Preserving the order is important.
	 */
	private static final String MenuTitles[] = {
		"New Main Frame", 
		"Duplicate Editor",
		"Workspace Hierarchy Browser",
		"Workspace File Browser",
		"View All Jobs",
		"View All Servers",
		"User Logs",
		"Programmer Logs",
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
		"MainFrame", 
		"DupTab",
		"DataSetView",
		"FileView",
		"Jobs",
		"Server",
		"UserLog",
		"ProgLog",
	};
	
	/**
	 * The strings for the action commands generated by the various menu items
	 * in the help menu. The list of values are organized in the same order as
	 * the ordinal values of the ActionType enumeration. Preserving the order 
	 * is important.
	 */
	private static final String ActionCmds[] = {
		"NewFrame", 
		"DupTab",
		"DATASET_TREE",
		"DATASET_FILE",
		"JOB_LIST",
		"SERVER_LIST",
		"USER_LOGS",
		"PROGRAMMER_LOGS",
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
		"Open a new main frame to provide extra screen real estate to view data",
		"Add a duplicate view of the currently selected view to the same main frame",
		"View files in workspace in a hierarchy based on data sets",
		"View files in current workspace in as a simple but comprehensive list",
		"View the current and completed jobs in the current workspace",
		"A tabular view of all the servers currently configured in this workspace",
		"Opens a tabular view of logs that provide additional information to a user",
		"View a textual form of the programmer logs with internal logs",
	};
}
