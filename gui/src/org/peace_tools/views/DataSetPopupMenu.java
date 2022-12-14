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
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.views;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import org.peace_tools.core.DeleteDialog;
import org.peace_tools.core.MainFrame;
import org.peace_tools.core.job.clustering.ClusteringJobWizard;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.FileEntry.FileEntryType;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobBase;
import org.peace_tools.workspace.Workspace;

/**
 * Class to generate and handle pop-up menus for data set viewers.
 * 
 * This is a common pop-up menu class that is shared by both the
 * DataSetTreeView and DataSetFileListView classes to generate
 * a pop-up menu when the user left-clicks on a entry. This 
 * class also handles some of the common menu options between both
 * views.
 */
public class DataSetPopupMenu extends JPopupMenu implements ActionListener {
	/**
	 * Creates all the menu entries to be displayed in the menu.
	 * 
	 * @param mainFrame The main frame that logically owns this menu.
	 * 
	 * @param treeView If this flag is true, then it is assumed that
	 * this popup menu will be used with a data set tree view. This
	 * determines some of the options that are added to the menu.
	 */
	public DataSetPopupMenu(MainFrame mainFrame, boolean treeView) {
		this.mainFrame  = mainFrame;
		this.isTreeView = treeView;
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Generate MST & Clustering", "mst", this, 
				"images/16x16/MST.png", null, false, false));
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Generate MST & Clustering", "mst", this, 
				"images/16x16/Cluster.png", null, true, false));
		addSeparator();
		//--------------------------------------------------
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Show with default viewer", "view", this, 
				"images/16x16/DefaultView.png", null, true, false));
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Clustering summary graph", "summary", this, 
				"images/16x16/ClusterSummary.png", null, true, false));
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Overlap analysis graph", "overlap", this, 
				"images/16x16/OverlapView.png", null, true, false));
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Show with text editor", "text", this, 
				"images/16x16/TextView.png", null, true, false));
		addSeparator();
		//--------------------------------------------------
		add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
				"Remove entry from Workspace", "delete", this, 
				"images/16x16/Delete.png", null, true, false));

		if (treeView) {
			add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
					"Collapse entries in sub-tree", "collapse", this, 
					null, null, true, false));
			add(Utilities.createMenuItem(Utilities.MenuItemKind.MENU_ITEM, 
					"Expand entries in sub-tree", "expand", this, 
					null, null, true, false));
		}
	}

	/**
	 * Enable/disable menu items based on item selection.
	 * 
	 * This method is a convenience method that is used to enable or
	 * disable menu items depending on the user's selection.
	 *  
	 * @param entry The Workspace object that has been selected.
	 * @param isExpanded This parameter is meanigful only for treeView
	 * objects. This parameter indicates if the selected object is
	 * in expanded or collapsed state.
	 */
	public void updateItems(Object entry, boolean isExpanded) {
		// Setup flags to indicate if we are dealing with a MST file
		// entry or a CLS file entry to streamline this method.
		// boolean isMST = false;
		boolean isCLS = false;
		String  jobID = null;
		if (entry instanceof FileEntry) {
			FileEntry fe = (FileEntry) entry;
			// isMST = fe.getType().equals(FileEntryType.MST);
			isCLS = fe.getType().equals(FileEntryType.CLS);
			jobID = fe.getGFL().getJobSummary().getJobID();
		}
		getComponent(1).setEnabled(entry instanceof DataSet);
		// getComponent(1).setEnabled(entry instanceof MSTData);
		
		// Enable/display default and text views 
		final boolean isGFL = (entry instanceof GeneratedFileList);
		getComponent(3).setEnabled(!isGFL);
		getComponent(6).setEnabled(!isGFL);
		
		getComponent(4).setEnabled(isCLS);
		getComponent(5).setEnabled(isCLS);
		// Enable/disable expand or collapse options
		if (isTreeView) {
			getComponent(9).setEnabled(isExpanded);
			getComponent(10).setEnabled(!isExpanded);
		}
		// Disable the delete option if the job associated with this entry
		// is not yet finished running.
		boolean jobDone = true;
		Job job = null;
		if (jobID != null) {
			job = Workspace.get().getJobList().getjob(jobID);
		}
		if (job != null) {
			// Check if the job is done running.
			jobDone = JobBase.JobStatusType.FAILED.equals(job.getStatus()) ||
				JobBase.JobStatusType.SUCCESS.equals(job.getStatus());
		}
		// Update deletion status based on job status.
		getComponent(8).setEnabled(jobDone);
		// Save current entry for use in actionPerformed method
		this.entry = entry;
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// Setup flags to indicate if we are dealing with a MST file
		// entry or a CLS file entry to streamline this method.
		boolean isCLS = false;
		if (entry instanceof FileEntry) {
			FileEntry fe = (FileEntry) entry;
			isCLS = fe.getType().equals(FileEntryType.CLS);
		}
		// Perform various operations based on action command
		String cmd = arg0.getActionCommand();
		if (("view".equals(cmd)) && (entry != null)) {
			 // Check and display the view (if one does not exist)
	        mainFrame.getViewFactory().createView(entry, false, false);
		} else if (("text".equals(cmd)) && (entry != null)) {
			 // Check and display a text view (if one does not exist)
	        mainFrame.getViewFactory().createView(entry, false, true);
		} else if ("summary".equals(cmd) && (entry != null) && (isCLS)) {
			// Check and display cluster summary 
	        mainFrame.getViewFactory().createSummaryView((FileEntry) entry);
		} else if ("overlap".equals(cmd) && (entry != null) && (isCLS)) {
			// Check and display overlap view
	        mainFrame.getViewFactory().createOverlapView((FileEntry) entry);
		} else if ("delete".equals(cmd) && (entry != null)) {
			// Delete the entry via the delete dialog.
			DeleteDialog delDiag = new DeleteDialog(mainFrame, entry);
			delDiag.setVisible(true);
		} else if ("mst".equals(cmd)) {
			 // Display the wizard to create a new data set
			 ClusteringJobWizard jobWizard = 
				 new ClusteringJobWizard("Create New MST", mainFrame);
			 Utilities.centerPanel(mainFrame, jobWizard);
			 jobWizard.showWizard("http://www.peace-tools.org/");
		}
	}

	/**
	 * Add a given action listener to all menu items.
	 * 
	 * This is a helper method that can be used to add a given
	 * action listener to all the menu items in this pop up menu.
	 * 
	 * @param listener The listener to be added to each menu item
	 * in the pop up menu.
	 */
	public void addActionListener(ActionListener listener) {
		final int compCount = getComponentCount();
		for(int i = 0; (i < compCount); i++) {
			Component c = this.getComponent(i);
			if (c instanceof JMenuItem) {
				JMenuItem menu = (JMenuItem) c;
				menu.addActionListener(listener);
			}
		}
	}
	
	/**
	 * Convenient reference to the main frame class that logically owns
	 * this menu. This value is set in the constructor
	 * and is never changed.
	 */
	private final MainFrame mainFrame;

	/**
	 * This instance variable tracks the currently selected
	 * work space entry. This value is set via call to the 
	 * updateItems method in this class.
	 */
	private Object entry;
	
	/**
	 * Instance variable that tracks if this pop-up menu is being
	 * used with a tree view or with a table view.
	 */
	private final boolean isTreeView;
	
	/**
	 * A generated serialization UID to keep compiler happy. 
	 */
	private static final long serialVersionUID = -4789513083096044919L;

}
