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

import javax.swing.Box;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JToolBar;

import org.peace_tools.core.dataset.DataSetWizard;
import org.peace_tools.core.job.JobWizard;
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
 * @note This class extends the WindowAdapter so that it can intercept and
 * process window closing event. This window closing event is generated 
 * either when the user clicks the close button (typically at the top-right
 * corner) or chooses the "Quit" option
 */
public class FileMenuHelper extends WindowAdapter implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "File" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public FileMenuHelper(MainFrame mainFrame) {
		this.mainFrame = mainFrame;
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
	 */
	public JMenu createFileMenu(JToolBar toolbar) {
		// Create the actual menu.
		JMenu fileMenu = new JMenu("File  ");
		// Create menus to save and switch work space.
		JMenuItem item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, 
					"Save current workspace", 
					"Save the current status and associated information about this workspace",
					"Save", this, "images/24x24/SaveWorkspace.png", 
					null, true, false);
		fileMenu.add(item);
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
					"Switch workspace", 
					"Save current workspace and switch to a new workspace",
					"Switch", this, "images/24x24/SwitchWorkspace.png", 
					null, false, false);
		fileMenu.add(item);
		fileMenu.addSeparator();
		
	    // First create and add the new file creation menu options.
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"New Data Set", 
				"Adds a new data set with a EST file to the current workspace",
				"NewDataSet", this, "images/24x24/NewDataSet.png", 
				null, true, false);
		fileMenu.add(item);
		// Option to compute new MST.
		item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Compute MST",
				"Starts wizard to compute MST and clustering for a ESTs in a data set",
				"NewMST", this, "images/24x24/NewMST.png", 
				null, false, false);
		fileMenu.add(item);
		// Option to compute new Clustering data
		item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Compute Clusters",
				"Starts wizard to compute MST and clustering for a ESTs in a data set",
				"NewCluster", this, "images/24x24/NewCluster.png", 
				null, true, false);
		fileMenu.add(item);
		fileMenu.addSeparator();
		// Options to import a data set from a different work space.
		item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Import Data Set",
					"Import a data set from another workspace into current workspace",
				"Import", this, "images/24x24/Import.png", 
				null, false, false);
		fileMenu.add(item);
		fileMenu.addSeparator();
		// Finally add the quit option.
		item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Exit PEACE",
				"Save workspace information and exit out of PEACE GUI",
				"Exit", this, "images/24x24/Exit.png", 
				null, true, false);
		fileMenu.add(item);
		
		if (toolbar != null) {
			// Add some of the tools that we anticipate users to work
			// with frequently to the toolbar. Ideally we could permit the
			// user to choose the tools that go on the toolbar. Maybe this
			// feature can be included if we have sufficient demand for it.
			toolbar.add(Utilities.createToolButton("images/24x24/SaveWorkspace.png", 
					null, "Save", this, 
					"Saves current workspace information to disk", true));
			toolbar.add(Box.createHorizontalStrut(5));
			toolbar.add(Utilities.createToolButton("images/24x24/NewDataSet.png", 
					null, "NewDataSet", this, 
					"Add a new data set (includes EST data and " +
					"cluster data) to this workspace", true));
			toolbar.add(Utilities.createToolButton("images/24x24/NewMST.png", 
					null, "NewMST", this, 
					"Compute new MST & Cluster data using a given FASTA file " +
					"from an existing data set", false));
			toolbar.add(Utilities.createToolButton("images/24x24/NewCluster.png", 
					null, "NewCluster", this, 
					"Compute new MST & Cluster data using a given FASTA file " +
					"from an existing data set", true));
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
			 dsWizard.showWizard("http://www.peace-tools.org/");
		 }
		 if ("NewMST".equals(cmd) || "NewCluster".equals(cmd)) {
			 // Display the wizard to create a new data set
			 JobWizard jobWizard = 
				 new JobWizard("Create New MST", mainFrame);
			 Utilities.centerPanel(mainFrame, jobWizard);
			 jobWizard.showWizard("http://www.peace-tools.org/");
		 }
		 if ("Exit".equals(cmd) && (reallyQuit())) {
			 System.exit(0);
		 }
		 if ("Save".equals(cmd)) {
			 mainFrame.saveDelayedWorkspace();
		 }
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
	 * Override the default method in WindowAdapter to process window closing.
	 * 
	 * This method is invoked because this class is added as a listener to
	 * the main frame (see constructor). When the user attempts to close 
	 * the main frame, this method intercepts this request and confirms with
	 * the user before closing the window.
	 * 
	 * @param e The window event associated with this call back. Currently this
	 * event is not used.
	 */
	@Override
	public void windowClosing(WindowEvent e) {
		if (reallyQuit()) {
			System.exit(0);
		}
	}
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
}
