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
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.AbstractButton;
import javax.swing.Box;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionListener;

import org.peace_tools.core.dataset.DataSetWizard;
import org.peace_tools.generic.Utilities;

/**
 * The file menu helper for the application.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "File" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation  for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that the FileMenuHelper is essentially an event
 * handler that is set on the various menu items in this class. This helper
 * class also provides a  createFileMenu() method that actually creates
 * the file menu.
 * 
 * <p><b>Note:</b>  This class adds a WindowAdapter so that it can intercept and
 * process window closing event. This window closing event is generated 
 * either when the user clicks the close button (typically at the top-right
 * corner) or chooses the "Quit" option.</p>
 */
public class FileMenuHelper extends AbstractMenuHelper implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "File" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public FileMenuHelper(MainFrame mainFrame) {
		super(HelperType.FILE_MENU, mainFrame);
		// Add a window closing handler to confirm exit.
		mainFrame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				if (reallyQuit()) {
					System.exit(0);
				}
			}
		});
	}
	
	/**
	 * Helper method to create the file menu.
	 * 
	 * This method performs the actual task of creating the "File" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "File" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "File" menu.
	 * 
	 * @param toolbar The toolbar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the toolbar is null, then
	 * shortcuts are not added.
	 * 
	 * @param jobHelper The job menu helper that must be used to add job
	 * creation entries to the new menu in the "File" menu. 
	 */
	public JMenu createFileMenu(JToolBar toolbar, AbstractMenuHelper jobHelper) {
		// Create the actual menu.
		JMenu fileMenu = new JMenu("File  ");
		// Add entries to the new sub menu.
		fileMenu.add(getMenuItem(ActionType.NEW_DATASET, true));
		fileMenu.add(jobHelper.getMenuItem(ActionType.COMPUTE_MST, true));
		fileMenu.add(jobHelper.getMenuItem(ActionType.COMPUTE_CLUSTERS, true));
		fileMenu.addSeparator();
		//----------------------------------------------------------
		
		// Create menus to save and switch work space.
		fileMenu.add(getMenuItem(ActionType.SAVE_WORKSPACE, true));
		fileMenu.add(getMenuItem(ActionType.IMPORT_DATASET, true));
		fileMenu.addSeparator();
		fileMenu.add(getMenuItem(ActionType.EXIT_PEACE, true));
		
		if (toolbar != null) {
			// Add some of the tools that we anticipate users to work
			// with frequently to the toolbar. Ideally we could permit the
			// user to choose the tools that go on the toolbar. Maybe this
			// feature can be included if we have sufficient demand for it.
			toolbar.add(getTool(ActionType.SAVE_WORKSPACE, true));
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(getTool(ActionType.NEW_DATASET, true));
		}
		return fileMenu;
	}

	@Override
	public void actionPerformed(ActionEvent event) {
		 String cmd = event.getActionCommand();
		 if ("NewDataSet".equals(cmd)) {
			 // Display the wizard to create a new data set
			 DataSetWizard dsWizard = 
				 new DataSetWizard("Create Data Set", mainFrame);
			 Utilities.centerPanel(mainFrame, dsWizard);
			 dsWizard.showWizard("http://www.peace-tools.org/downloads/manual.pdf#page=28");
		 } else if ("Exit".equals(cmd) && (reallyQuit())) {
			 System.exit(0);
		 } else if ("Save".equals(cmd)) {
			 mainFrame.saveDelayedWorkspace();
		 }
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
	public JMenuItem getMenuItem(ActionType actionType, boolean mainMenu) {
		int index = actionType.ordinal() - ActionType.SAVE_WORKSPACE.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainMenu ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";
		// Create and return the main menu item
		JMenuItem item = 
			Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, MenuTitles[index],
				(mainMenu ? MenuSubTitles[index] : null),
				ActionCmds[index], this, IconPath, 
				null, true, false);
		// Track context sensitive items.
		if (actionType.ordinal() >= ActionType.OPEN_DEFAULT_VIEW.ordinal()) {
			contextItemList.add(item);
			Utilities.setEnabled(item, false);
		}
		return item;
	}
	
	@Override
	public AbstractButton getTool(ActionType actionType, boolean mainToolBar) {
		int index = actionType.ordinal() - ActionType.SAVE_WORKSPACE.ordinal();
		if ((index < 0) || (index >= MenuTitles.length)) {
			// Unsupported option
			return null;
		}
		// Setup icon path depending on menu type
		final String IconPath = "images/" + (mainToolBar ? "24x24/" : "16x16/") 
			+ IconNames[index] + ".png";

		AbstractButton tool =
			Utilities.createToolButton(IconPath, null, ActionCmds[index], this, 
				MenuSubTitles[index], true);
		// Track context sensitive items.
		if (actionType.ordinal() >= ActionType.OPEN_DEFAULT_VIEW.ordinal()) {
			contextItemList.add(tool);
			tool.setEnabled(false);
		}
		return tool;
	}

	@Override
	public TreeSelectionListener getTreeSelectionListener(JTree tree) {
		// TODO Auto-generated method stub
		return null;
	}
	
    /** Double check before exiting.
     * 
     * This method is invoked when the user selects the "Quit" option
     * or clicks on the "X" button in the window.  This method pops
     * up a dialog box to double check the user really wants to quit.
     * 
     * @return If the user chooses to quit then this method returns 
     * true. Otherwise it returns false.
     */
	private boolean reallyQuit() {
		int result = 
			JOptionPane.showConfirmDialog(mainFrame, "Are you sure you want to exit?",
				"Double check", JOptionPane.YES_NO_OPTION);
		if (result == JOptionPane.YES_OPTION) {
			// Try and save the workspace before dropping out.
			mainFrame.saveWorkspace(null);
		}
		return result == JOptionPane.YES_OPTION;
	}
	
	/**
	 * The strings for each menu title created by this helper. The list of
	 * values are organized in the same order as the ordinal values of 
	 * the ActionType enumeration. Preserving the order is important.
	 */
	private static final String MenuTitles[] = {
		"Save current workspace", 
		"Switch workspace",
		"New Data Set",
		"Import Data Set",
		"Exit PEACE",
		"Open default viewer",
		"View cluster summary graph",
		"Open file in text viewer"
	};
	
	/**
	 * The icon file names for each menu item created by this helper. 
	 * The list of values are organized in the same order as the ordinal 
	 * values of the ActionType enumeration. Preserving the order is
	 * important. Note that a prefix directory (such as: images/16x16)
	 * and a suffix extension (.png) is added when tools or menu items
	 * are created.
	 */
	private static final String IconNames[] = {
		"SaveWorkspace", 
		"SwitchWorkspace",
		"NewDataSet",
		"Import",
		"Exit",
		"DefaultView",
		"ClusterSummary",
		"TextView"
	};
	
	/**
	 * The strings for the action commands generated by the various menu items
	 * in the help menu. The list of values are organized in the same order as
	 * the ordinal values of the ActionType enumeration. Preserving the order 
	 * is important.
	 */
	private static final String ActionCmds[] = {
		"Save", 
		"Switch",
		"NewDataSet",
		"Import",
		"Exit",
		"DefaultView",
		"ClusterSummary",
		"TextView"
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
		"Save the current status and associated information about this workspace", 
		"Save current workspace and switch to a new workspace",
		"Adds a new data set with a EST file to the current workspace",
		"Import a data set from another workspace into current workspace",
		"Save workspace information and exit out of PEACE GUI",
		"Open the selected file in its default viewer",
		"Show the selected clustering data file in a graphical view",
		"View the selected file in a plain text viewer"
	};
}
