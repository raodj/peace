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

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JToolBar;

import org.peace_tools.generic.Utilities;

/**
 * The file menu helper for the application.
 * 
 * This class encapsulates the code related to the operations performed by
 * various menu items in the "Job" menu. This class is typically created once
 * from the MainFrame.createMenus() method. The primary motivation  for
 * introducing sub-menu handler classes is to improve code organization and
 * minimize code clutter. Note that the JobMenuHelper is essentially an event
 * handler that is set on the various menu items in this class. This helper
 * class also provides a  createJobMenu() method that actually creates
 * the file menu.
 */
public class JobMenuHelper extends WindowAdapter implements ActionListener {
	/**
	 * The constructor. This class is an action listener that responds to the
	 * user clicking on various menu items. Since this is only an action listener
	 * the constructor does not have any special tasks to perform.
	 * 
	 * @param mainFrame The main frame that logically owns the "File" menu in its
	 * top-level menu bar. This reference is saved in this class for future use.
	 */
	public JobMenuHelper(MainFrame mainFrame) {
		this.mainFrame = mainFrame;
		assert(this.mainFrame != null);
	}
	
	/**
	 * Helper method to create the job menu.
	 * 
	 * This method performs the actual task of creating the "Job" menu. This
	 * method has been introduced to organize all the methods related to the
	 * "Job" menu into a single class. This method is invoked from the MainMenu
	 * class to create the "Job" menu.
	 * 
	 * @param toolbar The tool bar to which frequently used shortcuts can
	 * be added typically in the form of icons. If the tool bar is null, then
	 * shortcuts are not added.
	 * 
	 * @param fileHelper The action listener to use for handling items that
	 * are shared between file and job menu helper.
	 * 
	 * @param viewHelper The action listener to use for handling itesm that
	 * are shared between job and view menu.
	 */
	public JMenu createJobMenu(JToolBar toolbar, ActionListener fileHelper, 
			ActionListener viewHelper) {
		// Create the actual menu.
		JMenu jobMenu = new JMenu("Job  ");
		// Option to compute new MST.
		JMenuItem item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Job to Compute MST",
				"Starts wizard to schedule job to compute MST and clustering for a ESTs in a data set",
				"NewMST", fileHelper, "images/24x24/NewMST.png", 
				null, true, false);
		jobMenu.add(item);
		// Option to compute new Clustering data
		item = 
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Job to Compute Clusters",
				"Computes new clustering using existing MST data (quick operation)",
				"NewClusters", fileHelper, "images/24x24/NewCluster.png", 
				null, true, false);
		jobMenu.add(item);
		jobMenu.addSeparator();
		// Options to clear all jobs in the current work space.
		item =
			Utilities.createMenuItem(Utilities.MENU_ITEM, "Clear old jobs",
					"Clear all currently completed jobs in this workspace",
					"ClearJobs", this, "images/24x24/DeleteOldJobs.png", 
					null, true, false);
		jobMenu.add(item);
		jobMenu.addSeparator();
		//----------------------------------------------------------
		item = Utilities.createMenuItem(Utilities.MENU_ITEM, 
				"View All Jobs",
				"View the current and completed jobs in the current workspace",
				"ViewJobs", viewHelper, "images/24x24/Jobs.png", null, true, false);
		jobMenu.add(item);
		return jobMenu;
	}

	@Override
	public void actionPerformed(ActionEvent event) {

	}

	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu in its JMenuBar. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;
}
